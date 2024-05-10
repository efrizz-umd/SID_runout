% volumeFraction_updateFinderChunk_overlaps.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [volfracbins, Fmagbins, VYbins, Fxbins, Fybins, Fzbins, VXbins, VZbins, Vmagbins, xlabels, ylabels, bulkdensity,heights,...
    delta_bins, bins_stress_xx, bins_stress_yy, bins_stress_zz, bins_delta_total, ...
    bins_stress_xy, bins_stress_xz, bins_stress_yz,bins_particle_ids, ...
    store_psensor_ID, ...
    store_psensor_Fx, store_psensor_Fy, store_psensor_Fz, ...
    store_psensor_x, store_psensor_y, store_psensor_z, ...
    store_psensor_VY, store_psensor_VX, store_psensor_VZ, ...
    store_psensor_delta, store_psensor_heights, ...
    store_psensor_stress_xx, store_psensor_stress_yy, store_psensor_stress_zz, ...
    store_psensor_stress_xy, store_psensor_stress_xz, store_psensor_stress_yz, voro_bins,voro_bins_end] ...
               = volumeFraction_updateFinderChunk_overlaps(filename,filename_overlaps, ...
               spacing,ystart,ystop,zstart,zstop,startfilecheck,bin_particles_tocheck, ...
               store_psensor_ID, ...
               store_psensor_Fx, store_psensor_Fy, store_psensor_Fz, ...
               store_psensor_x, store_psensor_y, store_psensor_z, ...
               store_psensor_VY, store_psensor_VX, store_psensor_VZ, ...
               store_psensor_delta, store_psensor_heights, ...
               store_psensor_stress_xx, store_psensor_stress_yy, store_psensor_stress_zz, ...
               store_psensor_stress_xy, store_psensor_stress_xz, store_psensor_stress_yz,end_tstep_ind,SF)


% ************************************************************************
% This function opends the data file and reads in the data. The function the
% finds the binned (particle-averaged) data for a new pair of post and post_computes
% files (corresponding to a specific time step) and updates the storage arrays (appends the new data).
% The arrays are passed back and forth from this script and the wrapper script (there is definitely a
% more efficient way to do this).
%
% We describe generally the inputs and outputs instead of detailing them exahuastively
% they are all either either cell or double arrays
%
% inputs
% - filename - file name corresponding to particle states (ex: dump0.post)
% - filenamerun_overlaps - file name corresponding to overlaps (ex: dump_computes.post)
% - spacing - grid spacing (m)
% ystart, ystop, zstart, zstop - channel radial bounds and depth bounds
% - store_psensor_* - similar to bins (see outputs), but corresonds to the designated
% sensor particles (not presently used)
%
% file - file location
% step - the time step of the file to consider
% spacing - grid spacing to use (m)
% initflag - signals if this run is at the first time step, if so, set
% gridbounds off of this initial configuration
%
% output
% - *bins - is the storage array for a given quantity (volfrac = volume fraction, F = force
% V = velocity, stress = particle stresses, voro = bins of packing fraction computed
% via voronoi tesselation). There will either be an associated direction (x, y, z)
% or magnitude ('mag') flag for each named particle.
% - xlabels, ylabels - channel radial distance or depths

% ************************************************************************


%% Find bounds
delimiter = ' ';
startRow = 6;
endRow = 8;
formatSpec = '%f%f%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
fileID = fopen(filename,'r');

% import simulation bounds from the header
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, ...
'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow-1, 'ReturnOnError', false);
fclose(fileID);

% Allocate imported array to column variable names. hiBound is the side
% with the wall, and hiBound(2) is the wall-dependent boundary
loBound = dataArray{:, 1};
hiBound = dataArray{:, 2};

x_lo = loBound(1);
x_hi = hiBound(1);

clearvars delimiter startRow formatSpec fileID dataArray ans;

%% Read particle positions
% define file delimiters
delimiter = ' ';
startRow = 10;
% asterisk after the delimiter marker (%) means skip that field
formatSpec = '%f%*s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

% open file
fileID = fopen(filename,'r');

% read in data
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
id = dataArray{:, 1};
x1 = dataArray{:, 2};
y1 = dataArray{:, 3};
z1 = dataArray{:, 4};
vx = dataArray{:, 5};
vy = dataArray{:, 6};
vz = dataArray{:, 7};
r1 = dataArray{:, 8};
fx = dataArray{:, 9};
fy = dataArray{:, 10};
fz = dataArray{:, 11};
m = dataArray{:, 12};
stress_xx = dataArray{:, 13};
stress_yy = dataArray{:, 14};
stress_zz = dataArray{:, 15};
stress_xy = dataArray{:, 16};
stress_xz = dataArray{:, 17};
stress_yz = dataArray{:, 18};
% stress goes xx yy zz xy xz yz

% evaluation region
evallim = [x_lo,ystart,zstart;
            x_hi,ystop,zstop];


% compute overlaps
% average overlap is returned in an array, with array position
% corresponding to particle ID
[delta,delta_total] = overlap_Finder_PID(filename_overlaps,max(id));


% launch force update function here
[volfracbins, Fmagbins, Fxbins, Fybins, Fzbins, VXbins, VYbins, VZbins, Vmagbins, xlabels,ylabels,heights, ...
    delta_bins, bins_stress_xx, bins_stress_yy, bins_stress_zz, bins_delta_total, ...
    bins_stress_xy, bins_stress_xz, bins_stress_yz,bins_particle_ids, ...
    store_psensor_ID, ...
    store_psensor_Fx, store_psensor_Fy, store_psensor_Fz, ...
    store_psensor_x, store_psensor_y, store_psensor_z, ...
    store_psensor_VY, store_psensor_VX, store_psensor_VZ, ...
    store_psensor_delta, store_psensor_heights, ...
    store_psensor_stress_xx, store_psensor_stress_yy, store_psensor_stress_zz, ...
    store_psensor_stress_xy, store_psensor_stress_xz, store_psensor_stress_yz,voro_bins,voro_bins_end]= ...
         force_updateChunkFinder_overlaps(x1,y1,z1,r1,vx,vy,vz,fx,fy,fz,m,spacing,1,evallim,delta, delta_total, ...
         stress_xx, stress_yy, stress_zz, stress_xy, stress_xz, stress_yz,id,startfilecheck,bin_particles_tocheck, ...
         store_psensor_ID, ...
         store_psensor_Fx, store_psensor_Fy, store_psensor_Fz, ...
         store_psensor_x, store_psensor_y, store_psensor_z, ...
         store_psensor_VY, store_psensor_VX, store_psensor_VZ, ...
         store_psensor_delta, store_psensor_heights, ...
         store_psensor_stress_xx, store_psensor_stress_yy, store_psensor_stress_zz, ...
         store_psensor_stress_xy, store_psensor_stress_xz, store_psensor_stress_yz,end_tstep_ind,x_lo,x_hi,SF);


% find the bulk density
depth = .01; % ignore the top two centimeters
totalvolume = abs(evallim(2,1) - evallim(1,1))*abs(ystop - ystart)*abs((zstop-depth) - zstart);
particlemass = 0;
L = length(m);

% mass calculation
for jj = 1:L
    px = x1(jj);
    py = y1(jj);
    pz = z1(jj);
    % increment mass if within bounds
    if (px < evallim(2,1)) && (px > evallim(1,1)) && (py < ystop) && (py>ystart) ...
            && (pz < (zstop - depth)) && (pz > zstart)

        particlemass = particlemass + m(jj);
    end
end

bulkdensity = particlemass/totalvolume;
