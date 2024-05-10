% HarvestNoPlot_Production.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [grainstring, channelstring, porstring] ...
    = HarvestNoPlot_Production(filename,filenamerun_computes,timesteparray, filename_store,harvestpath,simtimestep,wspacing,SF)

% ************************************************************************
% This function loops through LIGGGHTS post files (states and computed overlaps)
% and performs a binning process, storing the results in a *.mat file. This function
% also calculates bed height, local packing fraction, and bulk density. The ability to
% track designated sensor particles is also included, but not currently in use.
%
% inputs
% filename - file name corresponding to particle states (ex: dump0.post)
% filenamerun_computes - file name corresponding to overlaps (ex: dump_computes.post)
% filename_store - file name to assing to file *.mat file
% harvestpath - directory containing code and output files
% simtimestep - simulation timestep (in seconds)
% wspacing - size of grid spacing
% SF - scaling factor
%
% outputs (strings corresponding to calculated:)
% grainstring - initial bed height (cm)
% channelstring - channel lenght (m)
% porstring - packing fraction (percent)
% ************************************************************************

% find the grain line
[grainline, channellength] = findGrainLineEvery(filename,0,SF)

% time and date information and strings
datetime.setDefaultFormats('default','yyyy-MM-dd')
currenttime = datetime;
formatOut = 'mmddyyyy';
timestring = datestr(currenttime,formatOut);
grainstring = num2str(round(100*grainline));
channelstring = num2str(channellength);


storestring = filename_store;

%% harvest data


% original file spacing - vestigial
windowrange = [0.0, (channellength)];
timestep1 = 1000;
timestep2 = 10000;
startFile = 0;
endFile = 225000;
startFile2 = 230000;%

L = length(timesteparray);

count = 1;

% storage for inner loop
timeplot = zeros(L,1);
% details for outer loop
ynumgrids = floor((windowrange(1,2) - windowrange(1,1)) / wspacing);
znumgrids = floor(grainline/wspacing);


% initialize data storage
ssdilationind_binary = zeros(znumgrids,ynumgrids);
ssdilationind_full = zeros(znumgrids,ynumgrids);
maxdtstep = zeros(ynumgrids,znumgrids);
store_vf_mean_pre = cell(znumgrids,ynumgrids);
store_vf_err_pre = cell(znumgrids,ynumgrids);
store_vf_mean_post = cell(znumgrids,ynumgrids);
store_vf_err_post = cell(znumgrids,ynumgrids);
store_timeplot = cell(znumgrids,ynumgrids);
store_vf_pdiff = cell(znumgrids,ynumgrids);
store_Fmag = cell(znumgrids,ynumgrids);
store_Fx = cell(znumgrids,ynumgrids);
store_Fy = cell(znumgrids,ynumgrids);
store_Fz = cell(znumgrids,ynumgrids);
Fmag_max = zeros(znumgrids,ynumgrids);
Fmag_maxtime = zeros(znumgrids,ynumgrids);
store_VY = cell(znumgrids,ynumgrids);
store_VX = cell(znumgrids,ynumgrids);
store_VZ = cell(znumgrids,ynumgrids);
store_heights = cell(1,L);
vy_max = zeros(znumgrids,ynumgrids);
store_Vmag = cell(znumgrids,ynumgrids);
vy_maxtime = zeros(znumgrids,ynumgrids);
store_delta = cell(znumgrids,ynumgrids);
store_delta_total = cell(znumgrids,ynumgrids);
store_stress_xx = cell(znumgrids,ynumgrids);
store_stress_yy = cell(znumgrids,ynumgrids);
store_stress_zz = cell(znumgrids,ynumgrids);
store_stress_xy = cell(znumgrids,ynumgrids);
store_stress_xz = cell(znumgrids,ynumgrids);
store_stress_yz = cell(znumgrids,ynumgrids);
store_voro_bins = cell(znumgrids,ynumgrids);

% particle sensors
% Passed back and forth to the functions to get the over time details
store_psensor_ID = zeros(znumgrids,1); % this one is just logged initially
store_psensor_Fx = cell(znumgrids,1);
store_psensor_Fy = cell(znumgrids,1);
store_psensor_Fz = cell(znumgrids,1);
store_psensor_x = cell(znumgrids,1);
store_psensor_y = cell(znumgrids,1);
store_psensor_z = cell(znumgrids,1);
store_psensor_VY = cell(znumgrids,1);
store_psensor_VX = cell(znumgrids,1);
store_psensor_VZ = cell(znumgrids,1);
store_psensor_delta = cell(znumgrids,1);
store_psensor_heights = cell(znumgrids,1);
store_psensor_stress_xx = cell(znumgrids,1);
store_psensor_stress_yy = cell(znumgrids,1);
store_psensor_stress_zz = cell(znumgrids,1);
store_psensor_stress_xy = cell(znumgrids,1);
store_psensor_stress_xz = cell(znumgrids,1);
store_psensor_stress_yz = cell(znumgrids,1);
store_psensor_mag_normal = cell(znumgrids,1);
store_psensor_mag_shear = cell(znumgrids,1);


ystart = windowrange(1);
ystop = windowrange(2);
zstart = grainline - znumgrids*wspacing;
zstop = grainline;

startfilecheck = 1;
% scheduler loops

% harvest data loop
% initialize particle bins
bin_particles_tocheck = 0;
end_step_ind = 0;

for i = 1:L
    timesteparray(i)
    a = strcat('dump',strcat(num2str(timesteparray(i)),'.post'));
    b = strcat('dump_computes',strcat(num2str(timesteparray(i)),'.post'));

    % set end step ind as necessary
    if i == L
        end_step_ind = 1;
    end

    % open the specified file and perform binning
    [vfracbins, Fmagbins, VYbins, Fxbins, Fybins, Fzbins, VXbins, VZbins, Vmagbins,xlabels,ylabels,bulkdensity,heights, ...
        delta_bins, bins_stress_xx, bins_stress_yy, bins_stress_zz, bins_delta_total, ...
        bins_stress_xy, bins_stress_xz, bins_stress_yz,bins_particle_ids, ...
        store_psensor_ID, ...
        store_psensor_Fx, store_psensor_Fy, store_psensor_Fz, ...
        store_psensor_x, store_psensor_y, store_psensor_z, ...
        store_psensor_VY, store_psensor_VX, store_psensor_VZ, ...
        store_psensor_delta, store_psensor_heights, ...
        store_psensor_stress_xx, store_psensor_stress_yy, store_psensor_stress_zz, ...
        store_psensor_stress_xy, store_psensor_stress_xz, store_psensor_stress_yz,...
	       voro_bins,voro_bins_end] ...
      = volumeFraction_updateFinderChunk_overlaps(strcat(filename,a),strcat(filenamerun_computes,b), ...
      wspacing,ystart,ystop,zstart,zstop,startfilecheck,bin_particles_tocheck, ...
      store_psensor_ID, ...
      store_psensor_Fx, store_psensor_Fy, store_psensor_Fz, ...
      store_psensor_x, store_psensor_y, store_psensor_z, ...
      store_psensor_VY, store_psensor_VX, store_psensor_VZ, ...
      store_psensor_delta, store_psensor_heights, ...
      store_psensor_stress_xx, store_psensor_stress_yy, store_psensor_stress_zz, ...
      store_psensor_stress_xy, store_psensor_stress_xz, store_psensor_stress_yz,end_step_ind,SF);


      % Update strings on first time step
    if timesteparray(i) == startFile
        porosity = mean(vfracbins(:));
        porstring = num2str(round(100*porosity));
        startfilecheck = startfilecheck + 1;
    end

    % update for next loop
    bin_particles_tocheck = bins_particle_ids;




    % this orients the grid storage as desired, with depth 0, r = 0 in
    % position row 1 x column 1
    vfracbins = flipud(vfracbins.');
    Fmagbins = flipud(Fmagbins.');
    Fxbins = flipud(Fxbins.');
    Fybins = flipud(Fybins.');
    Fzbins = flipud(Fzbins.');
    VYbins = flipud(VYbins.');
    VXbins = flipud(VXbins.');
    VZbins = flipud(VZbins.');
    Vmagbins = flipud(Vmagbins.');
    delta_bins = flipud(delta_bins.');
    bins_delta_total = flipud(bins_delta_total.');
    bins_stress_xx = flipud(bins_stress_xx.');
    bins_stress_yy = flipud(bins_stress_yy.');
    bins_stress_zz = flipud(bins_stress_zz.');
    bins_stress_xy = flipud(bins_stress_xy.');
    bins_stress_xz = flipud(bins_stress_xz.');
    bins_stress_yz = flipud(bins_stress_yz.');
    voro_bins = flipud(voro_bins.');
    voro_bins_end = flipud(voro_bins_end.');

    % store data and perform initial computes from the binned data

    % loop
    for j = 1:ynumgrids % vertical bins
        for k = 1:znumgrids % laters bins

          % initialize loop variables
            if timesteparray(i) == startFile
	           %[j, k]
            %[a, b] = size(vfracbins)

            vftemp = zeros(1,L);
            vftemp(:) = vfracbins(k,j);
            store_vf_mean_pre{k,j} = vftemp;
            end

            % accumulate values
            store_vf_mean_post{k,j} = [store_vf_mean_post{k,j}, vfracbins(k,j)];


            % append the next time step data
            store_vf_pdiff{k,j} = [store_vf_pdiff{k,j}, ...
                100 * (vfracbins(k,j) - store_vf_mean_pre{k,j}(1))/store_vf_mean_pre{k,j}(1)];

            store_Fmag{k,j} = [store_Fmag{k,j}, Fmagbins(k,j)];
            store_Fx{k,j} = [store_Fx{k,j}, Fxbins(k,j)];
            store_Fy{k,j} = [store_Fy{k,j}, Fybins(k,j)];
            store_Fz{k,j} = [store_Fz{k,j}, Fzbins(k,j)];
            store_VY{k,j} = [store_VY{k,j}, VYbins(k,j)];
            store_VX{k,j} = [store_VX{k,j}, VXbins(k,j)];
            store_VZ{k,j} = [store_VZ{k,j}, VZbins(k,j)];
            store_Vmag{k,j} = [store_Vmag{k,j}, Vmagbins(k,j)];

            store_delta{k,j} = [store_delta{k,j}, delta_bins(k,j)];
	          store_delta_total{k,j} = [store_delta_total{k,j}, bins_delta_total(k,j)];
            store_stress_xx{k,j} = [store_stress_xx{k,j}, bins_stress_xx(k,j)];
            store_stress_yy{k,j} = [store_stress_yy{k,j}, bins_stress_yy(k,j)];
            store_stress_zz{k,j} = [store_stress_zz{k,j}, bins_stress_zz(k,j)];
            store_stress_xy{k,j} = [store_stress_xy{k,j}, bins_stress_xy(k,j)];
            store_stress_xz{k,j} = [store_stress_xz{k,j}, bins_stress_xz(k,j)];
            store_stress_yz{k,j} = [store_stress_yz{k,j}, bins_stress_yz(k,j)];
	          store_voro_bins{k,j} = [store_voro_bins{k,j}, voro_bins(k,j)];
        end


    end


    % store vs time heights
    store_heights{1,i} = heights;
    % statistics

    % storage
    timeplot(count,1) = timesteparray(i)*simtimestep;
    Ltp = length(timeplot);

    count = count + 1;

end

% flip the particle sensors to be in the correct orientation (depth ind 1
% = surface)
store_psensor_ID = flipud(store_psensor_ID);
store_psensor_Fx = flipud(store_psensor_Fx);
store_psensor_Fy = flipud(store_psensor_Fy);
store_psensor_Fz = flipud(store_psensor_Fz);
store_psensor_x = flipud(store_psensor_x);
store_psensor_y = flipud(store_psensor_y);
store_psensor_z = flipud(store_psensor_z);
store_psensor_VY = flipud(store_psensor_VY);
store_psensor_VX = flipud(store_psensor_VX);
store_psensor_VZ = flipud(store_psensor_VZ);
store_psensor_delta = flipud(store_psensor_delta);
store_psensor_heights = flipud(store_psensor_heights);
store_psensor_stress_xx = flipud(store_psensor_stress_xx);
store_psensor_stress_yy = flipud(store_psensor_stress_yy);
store_psensor_stress_zz = flipud(store_psensor_stress_zz);
store_psensor_stress_xy = flipud(store_psensor_stress_xy);
store_psensor_stress_xz = flipud(store_psensor_stress_xz);
store_psensor_stress_yz = flipud(store_psensor_stress_yz);


% save the data for speeding things up next time
save([harvestpath storestring])
