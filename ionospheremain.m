% ionospheremain.m - Generate a top-of-ionosphere distribution based on mirror force alone (for now)
% and feed the distribution into the PADIC algorithm.  Output a cdf file containing the simulated
% satellite spectra at the altitude indicated in the variable "s_sat".  This is a MATLAB implementation
% of a depreciated (for now) C++ version.
%
% In the terminology used here:
%      ==== magnetosphere (or mag) ====  Top of mirroring region (source population)
%
%
%      ====   satellite (or sat)   ====  Where output plots are generated at
%
%      ====   ionosphere (or ion)  ====  Top layer of PADIC, bottom of mirroring region
%      ====       Layers 1-XX      ====  PADIC layers  \/
%      ====       Layers 1-XX      ====                \/
%
% Version: 1.0, 01 April 2022
%
% Tom Gade, School of Physics and Astronomy, University of Minnesota, 2022



%Script portion

ILAT = 72.0;
RADIUS_EARTH = 6.3712e6;          %meters
s_ion = getSatAlt(ILAT, 620e3);   %alt 620 km
s_sat = getSatAlt(ILAT, 4000e3);  %alt 4000 km
s_mag = 19881647.2473464;         %to match c++ code exactly
B_ion = getBFieldAtS(s_ion);
B_sat = getBFieldAtS(s_sat);
B_mag = getBFieldAtS(s_mag);

load('./data/ionospherecharacteristics.mat');

E_dist = 10.^(0.5:4.0/95:4.5);        %Energy, Pitch Angle bins used by the simulation (hi res)
PA_dist = 0.0025:0.005:179.9975;
E_sat = 10.^(0.5:4.0/47:4.5);         %Energy, Pitch Angle bins used by the simulated satellite (lo res)
PA_sat = 5:10:175;

printtime("start");
mag_list = generateInitLists(E_dist, PA_dist, 4.0/95, 0.005);  %here ion_list is initial ion upg
satbins = satellite(mag_list, B_mag, B_ion, B_sat, E_sat, PA_sat);

ion_dng_list = mirrorListHiToLo(mag_list,B_mag,B_ion);  %here, ion_list is dng at ion

upgbins = scattermain(ion_dng_list, s_layer, B_layer, h_layer, p, Z, B_sat, E_dist, PA_dist, E_sat, PA_sat); %return upg bins ONLY
outbins = upgbins + satbins;

savefold = write_data(PA_sat, E_sat, outbins);
% save(strcat('runs/',savefold,'/',savefold,'.mat'));

printtime("end");


%End script portion


%%% Functions
% Below are the functions that the script runs off of.  These somewhat closely mirror the CPP
% versions; however they are made to operate on vectors/matrices and are therefore executed in less
% lines of code than CPP, and work a little differently in some cases.  Functions are based around
% 2 data structures: lists and bins.  "lists" are lists of particles that store energy, pitch angle,
% dN or dE, bin minimum pitch angle, and bin max pitch angle - particles represent a histogram-style
% "bin" of particle pitches and energies.  "bins" contain the dN or dE flux of the bin without storing
% the pitch angle or energy corresponding to a bin (these are stored in separate arrays).

function dist_sat = satellite(mag_list, B_mag, B_ion, B_sat, E_bins, PA_bins)
    %Computes dE flux seen at satellite from magnetospheric (outside the loss cone)
    %and ionospheric sources (inside the loss cone)
    
    %Takes in a list defined as:
    % list(:,1) = particle energies
    % list(:,2) = particle pitches (degrees) - PA bin center
    % list(:,3) = particle dN weights
    % list(:,4) = PA bin min
    % list(:,5) = PA bin max
    
    list_dng_sat = [];
    list_upg_sat = [];
    
    %%%% Get distribution at satellite and ionosphere from initial distribution at simulation top
    Aratio_mag_sat = sqrt(B_mag / B_sat);
    
    PAsat = mirrorPA(mag_list(:,2), B_mag, B_sat); %PA at satellite
    PAion = mirrorPA(mag_list(:,2), B_mag, B_ion); %PA at ionosphere
    dngmask = (PAsat > 0);  %those that don't reflect somewhere after the satellite
    
    upgmask = dngmask & (PAion < 0); %those that reflect prior to ion, but after sat 
                                     %subset of dngmask
    
    list_dng_sat(:,1) = mag_list(dngmask,1); %energies of downgoing particles
    list_dng_sat(:,2) = PAsat(dngmask); %mirrored PA of downgoing particles
    list_dng_sat(:,3) = mag_list(dngmask,3) * Aratio_mag_sat; %dN flux of downgoing particles
    
    list_upg_sat(:,1) = mag_list(upgmask,1); %energies of upgoing particles
    list_upg_sat(:,2) = 180-PAsat(upgmask); %distribution at satellite {mirror PA of upgoing parti
    list_upg_sat(:,3) = mag_list(upgmask,3) * Aratio_mag_sat; %dN flux of upgoing particles
    %%%%
    
    
    % Compute bin boundaries at satellite - changed by mirror force
    list_dng_sat(:,4) = mirrorPA(mag_list(dngmask,4), B_mag, B_sat); %lower bin bound at sat
    list_dng_sat(:,5) = mirrorPA(mag_list(dngmask,5), B_mag, B_sat); %upper bin bound at sat
    list_dng_sat(list_dng_sat(:,5) == -1, 5) = 90 - 1.0e-10; %set bin bounds to 90 if part of the bin spills inside the loss cone
    %Future work: Amend flux depending on how much of the bin spills over 90 deg, send the rest up
    %This shouldn't be a big deal now - bin size is usually 0.005 deg
    
    %upgmask is a proper subset of dngmask, so we use dngmask to size "upgmask" appropriately - that is
    %to the number of rows that list_dng_sat has - then we index for the electrons that end up upg
    list_upg_sat(:,4) = 180 - list_dng_sat(upgmask(dngmask),5);
    list_upg_sat(:,5) = 180 - list_dng_sat(upgmask(dngmask),4);
    
    if (any(any(mag_list(dngmask,:) < 0)))
        error("satellite: something in mag_list is negative.  This is not valid.");
    end
    
    %Compute ratio in steradians of a pitch angle bin at the top of the sim to sters at satellite
    ster_ratio_dng = steradian_ratio(mag_list, B_mag, B_sat);
    ster_ratio_upg = ster_ratio_dng(upgmask(dngmask)); %are the same as before at a const |dPA| from 90
    
    %dN flux * steradian ratio * energy - dE flux
    list_dng_sat(:,3) = list_dng_sat(:,3) .* ster_ratio_dng .* list_dng_sat(:,1);
    list_upg_sat(:,3) = list_upg_sat(:,3) .* ster_ratio_upg .* list_upg_sat(:,1);
    
    if (any(any(list_dng_sat(:,[2,4,5]) > 90)) || any(any(list_upg_sat(:,[2,4,5]) < 90)))
        error("satellite: downgoing particles have pitch > 90 or upgoing particles have pitch < 90.");
    end
    
    if (any(list_dng_sat(:,3) < 0) || any(list_upg_sat(:,3) < 0) )
        error("satellite: downgoing or upgoing fluxes are negative.");
    end
    
    dist_sat = binning2D(E_bins, PA_bins, list_dng_sat) + binning2D(E_bins, PA_bins, list_upg_sat);
end

function ion_top_upg = scattermain(ion_list, s_level, B_level, h_level, p, Z, B_sat, E_dist, PA_dist, E_sat, PA_sat)
    % Implements the main PADIC algorithm - layer-based Rutherford scattering
    % Tracks the percent of the distribution in a bin that has scattered / reflected, adjusts for
    % gyroradius ratio (due to magnetic field), calls "upgoing_level" which is responsible for each
    % layer's scattering/reflecting dynamics.
    
    if (any(ion_list(:,2) >= 90.0))
        error('scattermain: list of particles that are downgoing at ionosphere has upgoing particles!');
    end
    
    percent_scattered = zeros(size(ion_list,1),1); %size(x,1) == #rows == #particles
    ion_top_upg = zeros(max(size(PA_dist)), max(size(E_dist)));
    
    %for diagnostic plotting
    %lists_upg_layer = cell(1,max(size(s_level))-1);
    %nums_refl = lists_upg_layer;
    %nums_scat = nums_refl;
    %lists_pct_scat_layer = nums_refl;
    %lists_pct_scat_total = nums_refl;
    %
    
    for (level = 1:max(size(s_level))-1)
        fprintf('Layer: %i/%i, ',level,size(s_level)-1);
        fprintf('s: %.0f, B: %.6e\n',s_level(level),B_level(level));
        printtime(strcat('Layer ',num2str(level)));
        
        Aratio_level_ion = sqrt(B_level(level) / B_level(1));
        %also steradian ratio from level to top
        
        %old_pct_sct = percent_scattered; %for diagnostics - prior to this layer's upg processes
        [upg_level, percent_scattered, ~, ~ ] = ...                     %num_refl, num_scat] = ...
            upgoing_level(ion_list, percent_scattered, E_dist, PA_dist, B_level, h_level, level, p, Z);
        
        upg_level(:,3) = upg_level(:,3) * Aratio_level_ion;
        %lists_upg_layer{level} = upg_level; %for diagnostics - get data with flux adjusted to ion top
        upg_level(:,2) = mirrorPA(upg_level(:,2), B_level(level), B_level(1));
        %neglecting mirroring bin bounds, because it's no longer useful
        
        %nums_refl{level} = num_refl; %for diagnostics
        %nums_scat{level} = num_scat;
        %lists_pct_scat_layer{level} = percent_scattered - old_pct_sct;
        %lists_pct_scat_total{level} = percent_scattered;
        
        ion_top_upg = ion_top_upg + binning2D(E_dist, PA_dist, upg_level);
    end
    
    ion_top_list = binsToParticleList(E_dist,PA_dist,ion_top_upg,PA_dist(2)-PA_dist(1));     %create a list
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %ion_top_list(:,2) = mirrorPA(ion_top_list(:,2),B_level(1),B_sat);        %mirror to satellite
    [ ster_ratio_ion_sat, ion_top_list(:,4:5), ion_top_list(:,2), ~ ] = ...
        steradian_ratio(ion_top_list, B_level(1), B_sat);
    ion_top_list(:,3) = ion_top_list(:,3) .* ster_ratio_ion_sat .* sqrt(B_level(1) / B_sat);
    ion_top_upg = binning2D(E_sat, PA_sat, ion_top_list);              %bin to appropriate size (match c++)
    ion_top_upg = ion_top_upg.*E_sat;
    
    %diagnostic_plots_scatter_layer(ion_list, B_level, 620000:-580000/(max(size(s_level))-1):40000, ...
    %                               lists_upg_layer, nums_refl, nums_scat,...
    %                               lists_pct_scat_layer, lists_pct_scat_total,...
    %                               E_bins, PA_bins); %for diagnostics
end


%
%
% Layer reflect / scatter handling
%
%
function [ upg_level, pct_scat, num_refl, num_scat ] = ...
        upgoing_level(ion_top_dng_list, pct_scat, E_bins, PA_bins, B_level, h_level, level, p, Z)
    
    Aratio_ion_level = sqrt(B_level(1) / B_level(level));
    
    PA_lvl_top = mirrorPA(ion_top_dng_list(:,2), B_level(1), B_level(level));   %PA at top of lvl
    PA_lvl_btm = mirrorPA(ion_top_dng_list(:,2), B_level(1), B_level(level+1)); %PA at btm of lvl
    
    %These are distinct categories - either a particle is reflected or it's scattered (if it arrives at this layer)
    %If the particle reflects within the layer, it doesn't scatter
    refl = (PA_lvl_top >= 0 & PA_lvl_btm < 0 & pct_scat < 1.0); %scatter and reflect mask arrays
    scat = (PA_lvl_top >= 0 & PA_lvl_btm >=0 & pct_scat < 1.0);
    
    lvl_scat = zeros(size(pct_scat));
    upg_level = zeros(max(size(pct_scat),5));
    
    %Case where particles are reflected - particles have already had steradians accounted for
    %according to top of ionosphere - they reflect and end up at the ionosphere top upgoing at
    %a pitch angle "180 - downgoing PA", so steradian ratio is the same (reflected around 90)
    upg_level(refl,:) = ion_top_dng_list(refl,:);
    upg_level(refl,[2 4 5]) = 180.0 - upg_level(refl,[2 5 4]); %reflect center, min and max bin PAs
    upg_level(refl,3) = (1-pct_scat(refl)) * Aratio_ion_level .* upg_level(refl,3); %no ster ratio, as described above
    pct_scat(refl) = 1.0;
    
    if (any(upg_level(refl,:) < 0.0))
        error('upgoing_level: upg_level has negative values!');
    end
    
    %Case where particles are scattered - need to account for steradian ratio here
    upg_level(scat,:) = ion_top_dng_list(scat,:);  %these array locations are used to contain what will scatter - not upgoing particles
    [ ster_ratio_scat, upg_level(scat,4:5), upg_level(scat,2), ~ ] = ...
        steradian_ratio(ion_top_dng_list(scat,:), B_level(1), B_level(level));
    if (any(any(upg_level(scat,2) ~= PA_lvl_top(scat)))) %should be bitwise identical - same operations
        error("upgoing_level: PAs not the same!!");
    end
    
    
    for (i=1:max(size(Z)))  %scatter off each species and sum
        lvl_scat(scat) = lvl_scat(scat) + scatterpercent(pct_scat(scat), Z(i), p(level,i), ...
                                                         h_level(level), upg_level(scat,1:2));
    end
    
    %Adjust scatter percent to prevent scattering over 100% of downgoing distribution
    lvl_scat(lvl_scat>1.0) = 1.0;  %max 100% scatters - note: not valid when multiple scattering events is implemented
    over1 = (pct_scat+lvl_scat) > 1.0; %mask that checks if more than 100% has scattered incl this layer
    
    lvl_scat(over1) = 1.0 - pct_scat(over1); %if so, set to whatever is left to scatter (will make pct_scat = 1 in line 255)
    %lvl_scat(lvl_scat<0 & abs(lvl_scat) < 1e-8) = 0.0; %if floating point math results in a very (very) small neg num, set to 0
    
    pct_scat = pct_scat + lvl_scat;
    pct_scat(over1) = 1.0;  %so % scattered isn't over 1
    
    if (any(lvl_scat < 0.0))
        error('upgoing_level: lvl_scat has negative values!');
    end
    if (any(pct_scat > 1.0))
        error('upgoing_level: pct_scat has values over 1!');
    end
    
    %pct scattered * A ratio * dN flux * ster ratio
    upg_level(scat,3) = lvl_scat(scat) .* Aratio_ion_level .* upg_level(scat,3) .* ster_ratio_scat;
    
    %Turn downgoing particles into backscatter
    sct_down_bin = binning2D(E_bins,PA_bins,upg_level(scat,:));
    backscat_list = downwardToBackscatter(E_bins,PA_bins,sct_down_bin);
    
    [ ster_ratio_bsupg, ~, ~, ~ ] = ...
        steradian_ratio(backscat_list, B_level(level), B_level(1));
    backscat_list(:,3) = backscat_list(:,3) .* ster_ratio_bsupg; %account for ster ratio up to ion top
    
    upg_level = [upg_level(~scat,:); backscat_list];
    upg_level = upg_level(upg_level(:,2) ~= 0.0,:);
    
    num_scat = max(size(backscat_list)); %for diagnostics
    num_refl = sum(refl);
    if (num_refl + num_scat ~= max(size(upg_level)))
        error('upgoing_level: sum of nums does not equal array size!');
    end
end

function backscat_list = downwardToBackscatter(E_bin_ctrs, PA_bin_ctrs, data)
    % Turns a "bins" data structure indicating the dE of particles that are scattering
    % by pitch angle and energy into a list of upgoing particles due to Evans' aggregated
    % primary/secondary scattering plots, with log fits determined from other computational work
    
    %Average dE flux over PA Bins
    PAsize = max(size(PA_bin_ctrs));
    Esize  = max(size(E_bin_ctrs));
    
    dPA = PA_bin_ctrs(2) - PA_bin_ctrs(1); %assumes constant size dPA between bins
    dNsum = sum(data, 1); %sum across pitch bins for a given energy
    
    %Calculate upward dN flux (backscatter) per E bin (no angle data - summed)
    dNback = johnd_flux(E_bin_ctrs, dNsum);
    
    backscat_list = zeros(PAsize*Esize,5);
    %Distribute BS dN flux isotropically over pitch bins
    for (pa = 1:max(size(PA_bin_ctrs)))
        if (PA_bin_ctrs(pa) < 90.0)
            start = (pa-1)*Esize + 1; %start and stop indices
            stop = pa*Esize;
            backscat_list(start:stop,1) = E_bin_ctrs;
            backscat_list(start:stop,2) = 180 - PA_bin_ctrs(pa); %upgoing
            backscat_list(start:stop,3) = dNback / (PAsize / 2);
            backscat_list(start:stop,4) = 180 - PA_bin_ctrs(pa) - dPA/2;
            backscat_list(start:stop,5) = 180 - PA_bin_ctrs(pa) + dPA/2;
        end
    end
    backscat_list = backscat_list(backscat_list(:,3) ~= 0.0,:);
end

function dNbackscat = johnd_flux(E, dN_E)
    P_LOGM = 0.505;       %Primary log linefit values
    P_LOGB = -5.16;
    S_LOGM_LT25 = -0.975; %Secondary log linefit values
    S_LOGB_LT25 = -1.47;
    S_LOGM_GT25 = -1.95;
    S_LOGB_GT25 = -0.11;
    
    lowE = E <= 25.0;
    
    f = @(E,logm,logb) 10.^( logm * log10(E) + logb ); %combination of log values
    f_slo = @(E) f(E,S_LOGM_LT25,S_LOGB_LT25);  %secondary low E
    f_shi = @(E) f(E,S_LOGM_GT25,S_LOGB_GT25);  %secondary hi E
    f_pri = @(E) f(E,P_LOGM,P_LOGB);            %primary all E
    
    dNbackscat = zeros(size(E));
    for (Einc = 1:max(size(E))) %iterate over incident E       dNavg(1:Einc)  also had .* E(1:Einc)
        dNbackscat(1:Einc) = dNbackscat(1:Einc) + dN_E(Einc) .* E(Einc) .*  ...
                           ( f_slo(E(1:Einc)) .* lowE(1:Einc) ...                %secd <= 25
                           + f_shi(E(1:Einc)) .* ~lowE(1:Einc)...                %secd >  25
                           + f_pri(E(1:Einc)/E(Einc)) .* (10000.0/E(Einc)) );    
                           %+ f_pri(E(1:Einc)/E(Einc)) .* (10000.0./E(1:Einc)) ); %prim
    end
end


%
%
% Distribution generation and binning
%
%
function [mag_list, ion_list] = generateInitLists(E_dist, PA_dist, dlogE, dPA)
    % Generate initial lists of particles for mag and ion.  Consists of output matrices
    % of binsToParticleLists()
    % list(:,1) = particle energies
    % list(:,2) = particle pitches (degrees) - PA bin center
    % list(:,3) = particle dN weights
    % list(:,4) = PA bin min
    % list(:,5) = PA bin max

    % 
    % Maxwellian distribution with no variable kappa
    %
%     f_maxw = @(a1,a2,E_char_lo,E_char_hi,E) a1*exp(-E/E_char_lo) + a2*exp(-E/E_char_hi); %maxwellian distribution kappa = standard
%     f_mag_integ = @(E) 1./E .* f_maxw(4.5e4,157.5,4,4000, E); %integrating over bin gives total dN flux in bin
%     mag_maxw = zeros(1,max(size(E_dist)));
    
    %
    % Maxwellian distribution with variable kappa
    %
    kappa_func = @(a, kT, kappa, E) (a/kT*2/sqrt(pi).*E.*gamma(kappa+1)/kappa^1.5/gamma(kappa-.5).*E/kT/((1+E/(kT*kappa)).^(kappa+1)))./E./E;
    f_maxw = @(a1,a2,E_char_lo,E_char_hi, kappa, E) kappa_func(a1,E_char_lo, kappa, E) + kappa_func(a2, E_char_hi, kappa, E);  %maxwellian distribution variable kappa
    f_mag_integ = @(E) 1./E .* f_maxw(8e5,3e9,4,5000, 3, E); %integrating over bin gives total dN flux in bin
    mag_maxw = zeros(1,max(size(E_dist)));
    
    logE = log10(E_dist); % array of log10 of E_dist
    E_log_range = logE(end) - logE(1); % finds the range of the logE array (4)
    dlogE = E_log_range/max(size(E_dist)); % difference between distance of each E value in E_dist (.0417)
   
    
    E_bin_start = 10.^(logE - dlogE/2);
    
    E_bin_start = [E_bin_start.'; 10^(logE(end) + dlogE/2)];
    for (i = 1:max(size(E_dist))) % integrate 1/E * f_maxwellian(E) over bin range of E to get total dN flux
        mag_maxw(i) = integral(f_mag_integ, E_bin_start(i), E_bin_start(i+1));
    end

    disp(mag_maxw)
    
    scale = 6.1; %adjust fluxes to match real data
    

    %%%% Input Real MEPT Data
    fileID1 = fopen('./data/bins/particles_final/elec_vpara.bin');
    fileID2 = fopen('./data/bins/particles_final/elec_vperp.bin');
    elecParaData = fread(fileID1);
    elecPerpData = fread(fileID2);
    fclose(fileID1);
    fclose(fileID2);
    
    vel_mag = sqrt((elecParaData.^2)+(elecPerpData.^2));
    ergVal = .5.*(9.109e-28).*vel_mag.^2;
    eVVal = 6.242e26 * ergVal;
    
    Erange = [0,E_dist];
    real_mag_maxw = zeros(1,96);
    for i = 1:length(eVVal)
        j = 1;
        while true
            if eVVal(i) >= Erange(j) && eVVal(i) < Erange(j+1)
                real_mag_maxw(j) = real_mag_maxw(j)+1;
                break
            elseif eVVal(i) > Erange(97)
                real_mag_maxw(96) = real_mag_maxw(j)+1;
                break
            end
            j = j+1;
            
        end
    end
    
%     mag_maxw_real = readmatrix('./data/magnetospheremaxwellian.csv');  %read data from disk of maxwellian
%     mag_maxw_real = mag_maxw_real * scale;  %values adjusted to match real world fluxes
     
    if(false) %diagnostic plot - compare maxwellians with real world data
    figure; plot(E_dist, mag_maxw);
    set(gca,'XScale','log');
    set(gca,'YScale','log');
    hold on; plot(E_dist, mag_maxw_real);
    end
    
    if (true)  %use real distribution
%     mag_maxw = mag_maxw_real;
    mag_maxw = real_mag_maxw;
    end
    
    mask_dng = PA_dist <= 90.0; % creates 1x36000 vector of 1s if PA_dist is less than 90 and 0 otherwise
    size_dng = sum(mask_dng); % 18000
    mag_maxw = reshape(mag_maxw, [1,max(size(mag_maxw))]); % ensure shaped the right way
    mag_maxw = repmat(mag_maxw,size_dng,1); % reshapes 96x1 matrix to 18000x96 matrix where every column only has one value
    mag_list = binsToParticleList(E_dist,PA_dist(mask_dng),mag_maxw,dPA);
    
    %generate ionospheric list (if applicable)
    %size_upg = max(size(PA_dist)) - size_dng;
    %ion_maxw = readmatrix('./data/ionospheremaxwellian.csv');     %data derived from real sat data
    %ion_maxw = ion_maxw * 5.5;
    %ion_maxw = (1e3*exp(-E_dist/6));
    %ion_maxw = reshape(ion_maxw, [1,max(size(ion_maxw))]);
    
    %ion_maxw = repmat(ion_maxw,size_upg,1);
    
    %ion_list = binsToParticleList(E_dist,PA_dist(~mask_dng),ion_maxw);
    ion_list = [0,0,0]; %no ionospheric source particles - computed by using backscatter physics
end

function dist_ion = mirrorListHiToLo(init_list, B_fr, B_to)
    % Mirrors a list of particles from high altitude to low altitude
    % Discards any that reflect and returns a smaller sized "list" data structure
    
    %Takes in a list defined as:
    % list(:,1) = particle energies
    % list(:,2) = particle pitches (degrees) - PA bin center
    % list(:,3) = particle dN weights
    % list(:,4) = PA bin min
    % list(:,5) = PA bin max
    
    Aratio = sqrt(B_fr / B_to);
    
    [ ster_ratio, bounds, PAmir, mask ] = steradian_ratio(init_list, B_fr, B_to);
    
    dist_ion = [];
    dist_ion(:,1) = init_list(mask,1);
    dist_ion(:,2) = PAmir(mask);
    dist_ion(:,3) = init_list(mask,3) .* ster_ratio * Aratio;
    dist_ion(:,4) = bounds(mask,1);
    dist_ion(:,5) = bounds(mask,2);
    
    %dist_ion_old = [];
    %dist_ion_old(:,1) = init_list(:,1);
    %dist_ion_old(:,2) = mirrorPA(init_list(:,2), B_hi, B_lo);
    %dist_ion_old(dist_ion_old(:,2)>0,3) = init_list(dist_ion_old(:,2)>0,3) * Aratio;
end

function partlist = binsToParticleList(E_bins, PA_bins, binned_data, dPA)
    % Converts from "bins" data structure which contains dN or dE in each bin:
    % | dN(E1, PA1) or dE(E1, PA1) | dN(E2, PA1) or dE(E2, PA1) | dN(E3, PA1) or dE(E3, PA1) ....
    % | dN(E1, PA2) or dE(E1, PA2) | dN(E2, PA2) or dE(E2, PA2) | dN(E3, PA2) or dE(E3, PA2) ....
    % ......
    % for E = (energy low -> energy hi), evenly spaced in log space
    %     PA= (PA low -> PA hi), evenly spaced in decimal space
    %
    % To a "list" data structure:
    % list(:,1) = particle energies
    % list(:,2) = particle pitches (degrees) - PA bin center
    % list(:,3) = particle dN weights
    % list(:,4) = PA bin min
    % list(:,5) = PA bin max
    
    E_size = max(size(E_bins)); % 96
    PA_size = max(size(PA_bins)); % 18000
    
    if (any(size(binned_data) ~= [PA_size,E_size]))
        error("binsToParticleList: binned_data is malformed");
    end
    
    partlist = zeros(E_size*PA_size,5);
    for (i = 1:PA_size)
        for (j = 1:E_size)
            if (binned_data(i,j) ~= 0)
                partlist((i-1)*E_size+j,1) = E_bins(j);
                partlist((i-1)*E_size+j,2) = PA_bins(i);
                partlist((i-1)*E_size+j,3) = binned_data(i,j);
            end
        end
    end
    partlist = partlist(partlist(:,3) > 0, :);
    
    partlist(:,4) = partlist(:,2) - dPA/2;
    partlist(:,5) = partlist(:,2) + dPA/2;
end

function binned = binning2D(E_bin_ctrs, PA_bin_ctrs, list)
    % Turns a "list" data structure into a "bins" data structure
    %
    % REQUIRES EVENLY SPACED BINS (E in log space, PA in decimal space)
    
    % list(:,1) = particle energies
    % list(:,2) = particle pitches (degrees) - PA bin center
    % list(:,3) = particle dN weights
    % list(:,4) = PA bin min
    % list(:,5) = PA bin max
    
    if (size(list,2) ~= 5)
        error("binning2D: error data is malformed");
    end
    
    binned = zeros(max(size(PA_bin_ctrs)), max(size(E_bin_ctrs)));
    
    if (isempty(list)) %null list, zero flux in every bin
        return;
    end
    
    delta_E = log10(E_bin_ctrs(2)) - log10(E_bin_ctrs(1));
    delta_PA = PA_bin_ctrs(2) - PA_bin_ctrs(1);
    
    if (delta_E < 0 || delta_PA < 0)
        error('binning2D: E or PA bins are decreasing.  This isnt supported at the moment.');
    end
    
    E_bin_st = 10.^(log10(E_bin_ctrs) - 0.5*delta_E);
    E_bin_st(end+1) = 10^(log10(E_bin_st(end))+delta_E);
    PA_bin_st = PA_bin_ctrs - 0.5*delta_PA;
    PA_bin_st(end+1) = PA_bin_st(end)+delta_PA;
    
    %Get bin indices - this is the part that assumes evenly spaced bins
    list(:,6) = floor((log10(list(:,1))-log10(E_bin_st(1)))/delta_E) + 1;
    list(:,7) = floor((list(:,2)-PA_bin_st(1))/delta_PA) + 1;
    
    if (any(list(:,[6 7]) < 0))
        error('binning2D: bins with negative values were produced.  These are invalid!');
    end
    
    outcnt = 0;
    for (i = 1:size(list,1))
        if ( (list(i,6) > max(size(E_bin_ctrs))) || list(i,6) < 1 ||...
             (list(i,7) > max(size(PA_bin_ctrs)))|| list(i,7) < 1 )
            outcnt = outcnt + 1;
            continue;
        end
        if ( list(i,1) < E_bin_st(list(i,6))  || list(i,1) > E_bin_st(list(i,6)+1) || ...
             list(i,2) < PA_bin_st(list(i,7)) || list(i,2) > PA_bin_st(list(i,7)+1) )
            error("outside bin");
        end
        binned(list(i,7),list(i,6)) = binned(list(i,7),list(i,6)) + list(i,3);
    end
end


%
%
% Scattering / mirroring helper functions
%
%
function PAf_deg = mirrorPA(PA0_deg, Binit, Bfinal)
    % Mirrors an initial pitch angle due to the B field at Binit to Bfinal (due to mirror force only)
    % Produces similar values to those produced by MEPT - time-stepping particle dynamics code
    
    if ((any(PA0_deg < 90.0) && abs(Binit) > abs(Bfinal)) ||...
        (any(PA0_deg > 90.0) && abs(Binit) < abs(Bfinal)))
        error("mirrorPA: particle, Bfield gradient mismatch");
    end
    if (Binit == Bfinal)
        PAf_deg = PA0_deg;
        return;
    end
    
    %intermediate value for computing final angle
    PAf_deg = zeros(size(PA0_deg));
    PAf_deg(PA0_deg >=0) = Binit / Bfinal * (1 + 1./(tan_DEG(PA0_deg(PA0_deg >=0)).^2)) - 1;
    PAf_deg(PA0_deg < 0) = PA0_deg(PA0_deg < 0);
    
    PAf_deg(PAf_deg < 0) = -1; %if this is the case, particle reflects before Bfinal
    PAf_deg(PAf_deg >=0) = atan_DEG(sqrt(1./PAf_deg(PAf_deg >=0)));
    if (abs(Binit) > abs(Bfinal))  %if B field strength drops off, no mirroring occurs
        PAf_deg = 180 - PAf_deg;   %PA > 90 is defined as upgoing (away from Earth)
    end
end

function sctpct = scatterpercent(pctSctAbove, Z, p, h, partlist)
    % Computes the percentage of incident particles that scatter due to Rutherford scattering
    % For this initial model, if a particle scatters at all (even a little), it is counted in with
    % Evans-style aggregated primary / secondary curves, with log fits that deviate from Evans'
    % original values determined from other computational work
    
    if (max(size(partlist(:,1))) ~= max(size(pctSctAbove)))
        error("scatterpercentage: arrays are ill-formed");
    end
    E = partlist(:,1);
    PA = partlist(:,2);
    
    sctpct = (1-pctSctAbove) .* 1.62e-14 * Z^2 * p * h ./ (E.^2 .* cos_DEG(PA));
    
    sctpct(pctSctAbove >= 1.0) = 0;
    sctpct(PA > 90) = 0;
    sctpct(sctpct > 1) = 1;  %Not valid once partial scattering is considered
end

function [ ster_ratio, new_bounds, ctr_bin, mask ] = steradian_ratio(initlist, B_init, B_final)
    %Compute ratio in steradians of a pitch angle bin from some initial point to a final point
    %Since pitch angles are mirrored here, returns these to save other sections of code from recomputing
    
    new_bounds = zeros(size(initlist,1),2);
    new_bounds(:,1) = mirrorPA(initlist(:,4), B_init, B_final); %lower bin bound at final B
    new_bounds(:,2) = mirrorPA(initlist(:,5), B_init, B_final); %upper bin bound at final B

    ctr_bin = mirrorPA(initlist(:,2), B_init, B_final);
    mask = ctr_bin >= 0; %returns logical 1 if there is a nonreflected value
    if (any(new_bounds(mask,1) == -1))  %if this happens, this indicates an invalid value was passed in
        error("steradian_ratio: lower bin bound reflects");
    end
    new_bounds(new_bounds(mask,2) == -1, 2) = 90.0; %upper bound of the bin reflects - set upper bnd to 90
    
    ster_ratio = ...
    (cos_DEG(new_bounds(mask,1)) - cos_DEG(new_bounds(mask,2))) ./ (cos_DEG(initlist(mask,4))   - cos_DEG(initlist(mask,5)));
    %Future work: Amend flux depending on how much of the bin spills over 90 deg, send the rest up
    %This shouldn't be a big deal now - bin size is usually 0.005 deg
    
    %\/ This indicates a bin size of 0 before or after mirroring, which is invalid
    if (any(ster_ratio == 0) || any(isnan(ster_ratio)) || any(isinf(ster_ratio)))
        error("satellite: ratio of steradians is 0, NaN or inf.");
    end
    
end

function printtime(label)
    start = clock();
    start = fix(start);
    fprintf(label);
    fprintf(" : Time: %i:%i:%i\n\n", start(4), start(5), start(6));
end

function ret = sin_DEG(deg)
    RADS_PER_DEG = pi / 180;
    ret = sin(deg * RADS_PER_DEG); 
end

function ret = cos_DEG(deg)
    RADS_PER_DEG = pi / 180;
    ret = cos(deg * RADS_PER_DEG);
end

function ret = tan_DEG(deg)
    RADS_PER_DEG = pi / 180;
    ret = tan(deg * RADS_PER_DEG);
end

function deg = atan_DEG(num)
    RADS_PER_DEG = pi / 180;
    deg = atan(num) / RADS_PER_DEG;
end

function fname = write_data(PA_sat, E_sat, outbins)
    %Create CDF file with data
    fname = datestr(datetime(),'yymmdd.HHMM');
    mkdir(strcat('runs/',fname));
    cdfhandle = cdflib.create(strcat('runs/',fname,'/',fname,'.cdf'));

    szP = max(size(PA_sat));
    szE = max(size(E_sat));

    zE =  cdflib.createVar(cdfhandle, 'Mid-Bin Energies (eV)', 'CDF_DOUBLE', 1, szE,...
                           true, true);
    zA =  cdflib.createVar(cdfhandle, 'Mid-Bin Angles (Degrees)', 'CDF_DOUBLE', 1, szP,...
                           true, true);
    zdE = cdflib.createVar(cdfhandle, 'Electrons Energy/Pitch Angle Count', 'CDF_DOUBLE', 1,...
                           [ szE, szP ], true, [ true, true ]);

    cdflib.hyperPutVarData(cdfhandle, zE,  [0 1 1], {0 szE 1},   E_sat);
    cdflib.hyperPutVarData(cdfhandle, zA,  [0 1 1], {0 szP 1},   PA_sat);
    cdflib.hyperPutVarData(cdfhandle, zdE, [0 1 1], {[0 0] [szE szP] [1 1]}, outbins.');

    cdflib.close(cdfhandle);
    
    copyfile('ionospheremain.m',strcat('runs/',fname,'/ionospheremain.m'));
end