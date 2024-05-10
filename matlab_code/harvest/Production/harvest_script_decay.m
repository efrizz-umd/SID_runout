% harvest_script_decay.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

% This is a wrapper script that manages the data binning and accumulation process
% of LIGGGHTS output files. Files should be located in two separate directors located in
% the same folder as the wrapper script and utility functions (a 'dump' directory which
% contains the particles states and the 'dump_computes' which contains every particle-particle overlap).
% The output data should be set to the same rates as the code is written to function based off
% the assumption that there is a states and overlaps file for each timestep.

% * variables to be MANUALLY ADJUSTED *
% timestep1 - output rate, resolution 1 (ex: timestep1 = 500, if files output every 500 steps)
% timestep2 - output rate, resolution 2
% endFile1 - simulation step corresponding to the end of resolution 1 range
% endFile2 - simulation step corresponding to the end of the simulation
% simtimestep - 1 simulation timestep corresonds to this number of seconds
% SF - a scaling factor, corresonding to your monodisperse particle size relative to 1.25E-3 m
% filename_store - meaningful name that gives details about the test from which the binned data originated

% * results *
% Running the wrapper scripts launches HarvestNoPlot_Production.m, which generates a *.mat
% file with binned, vs-time data corresponding to particle velocites, forces, and stresses, as
% well as measuring bed height and bulk volume


% ------ the code ------ %

% for saving (filenames/filepaths)
datetime.setDefaultFormats('default','yyyy-MM-dd')
currenttime = datetime;
formatOut = 'mmddyyyy';
timestring = datestr(currenttime,formatOut);

thisdir = pwd;
filename = [thisdir '/post/'];
filename_computes = [thisdir '/post_computes/'];

% saved *.mat file name, MANUALLY ADJUSTED
filename_store = ['2m_channel_20cmfill_compact_10ms_shock_type_val_' timestring '.mat'];
workingdir = thisdir;
harvestpath = [workingdir '/mat_files/'];

filenamerun = convertStringsToChars(filename)
filenamerun_computes = convertStringsToChars(filename_computes)

% timestep breakdown, MANUALLY ADJUSTED
timestep1 = 250;
timestep2 = 2000;
startFile = 0;
endFile1 = 12000;
startFile2 = 14000;%
endFile2 = 500000;


% scaling factor, MANUALLY ADJUSTED
SF = 1000; % scaling factor corresponding to particle size table

% timestep array, MAUALLY ADJUSTED (if more than two output rates)
timesteparray = [startFile:timestep1:endFile1, startFile2:timestep2:endFile2];

% simulation time step, MANUALLY ADJUSTED
simtimestep = 0.00001; % time step within the LIGGGHTS inputscript

% grid spacing
wspacing = SF* 0.005; % grid spacig calculations

%% run the analysis function
[grainstring, channelstring, porstring] ...
    = HarvestNoPlot_Production(filenamerun,filenamerun_computes,timesteparray,filename_store,harvestpath,simtimestep,wspacing,SF);
