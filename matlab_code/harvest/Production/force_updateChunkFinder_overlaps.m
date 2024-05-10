% force_updateChunkFinder_overlaps.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [volfracbins, Fmagbins, Fxbins, Fybins, Fzbins, VXbins, VYbins, VZbins, Vmagbins, xlabelout, ylabelout,heights, ...
    delta_bins, bins_stress_xx, bins_stress_yy, bins_stress_zz, bins_delta_total, ...
    bins_stress_xy, bins_stress_xz, bins_stress_yz,bins_particle_ids, ...
    store_psensor_ID, ...
    store_psensor_Fx, store_psensor_Fy, store_psensor_Fz, ...
    store_psensor_x, store_psensor_y, store_psensor_z, ...
    store_psensor_VY, store_psensor_VX, store_psensor_VZ, ...
    store_psensor_delta, store_psensor_heights, ...
    store_psensor_stress_xx, store_psensor_stress_yy, store_psensor_stress_zz, ...
    store_psensor_stress_xy, store_psensor_stress_xz, store_psensor_stress_yz, ...
    voro_bins,voro_bins_end] ...
    = force_updateChunkFinder_overlaps(x1,y1,z1,r1,vx,vy,vz,fx,fy,fz,m,spacing,initflag, ...
    evalreg,overlaps, overlaps_total,stress_xx, stress_yy, stress_zz, stress_xy, stress_xz, stress_yz,p_id,startfilecheck,bin_particles_tocheck, ...
    store_psensor_ID, ...
    store_psensor_Fx, store_psensor_Fy, store_psensor_Fz, ...
    store_psensor_x, store_psensor_y, store_psensor_z, ...
    store_psensor_VY, store_psensor_VX, store_psensor_VZ, ...
    store_psensor_delta, store_psensor_heights, ...
    store_psensor_stress_xx, store_psensor_stress_yy, store_psensor_stress_zz, ...
    store_psensor_stress_xy, store_psensor_stress_xz, store_psensor_stress_yz,end_tstep_ind,x_lo,x_hi,SF)


% ************************************************************************
% This function establishes the grid bounds and manages updating the sensor particles
% as well as the running of the binning sripts. Also computes the channel height
%
% inputs
% x1, y1, z1 - particle positions (m)
% r1 - particle radius (m)
% vx, vy, vz - particle velocities (m/s)
% fx, fy, fz - particle forces (N)
% evalreg - grid boundaries
% overlaps, overlaps_total - arrays of the average and total overlap per particle
% stress_** - stress in a given direction (LIGGGHTS stress units)
% p_id - particle IDs
% filename - file name corresponding to overlaps (ex: dump_computes0.post)
% num particles - total number of particle in the sim
%
% outputs (strings corresponding to calculated:)
% - store_psensor_* - is the particle sensor storage array for a given quantity (volfrac = volume fraction, F = force
% V = velocity, stress = particle stresses, voro = bins of packing fraction computed
% via voronoi tesselation). There will either be an associated direction (x, y, z)
% or magnitude ('mag') flag for each named particle.
% ************************************************************************



pauseflag = 0;

% define end wall
wall = max(y1);

% set start and stop regions, can be modified later to be passed
ystart = .0; %
ystop = .0; %

% find the the particle line ('height')
zmax = max(z1) + 1.0001*r1(1);


% set evaluation region based on initial particle configuration if init
if initflag == 0

    % limits of evaluation region: [xlo, ylo, zlo; xhi, yhi, zhi]
    evallim = [x_lo,ystart,0;
                x_hi,(wall-ystop),zmax];
else
    evallim = evalreg;
end

% determine the grid y and z grid bounds. We average over x so there is only
% 1 x index
gBx = [evallim(1,1), evallim(2,1)];
gBy = (evallim(1,2):spacing:evallim(2,2));
gBz = (evallim(1,3):spacing:evallim(2,3));

% reset the end value, this will stretch the last grid a bit, but will make
% sure no area is left off
gBy(end) = evallim(2,2);
gBz(end) = evallim(2,3);


xlabelout = mean([gBy(1:end-1);gBy(2:end)]);
ylabelout = mean([gBz(1:end-1);gBz(2:end)]);

numgridy = length(gBy) - 1;
numgridz = length(gBz) - 1;

gridvolbins = zeros(numgridy,numgridz);

% sensor indicator to use, half way through channel
s_i = round(numgridy/2);

for j = 1:numgridy
    for k = 1:numgridz
        % find cube volume based on grid bounds
        gridvolbins(j,k)=abs(gBx(1) - gBx(2)) * ...
                        abs(gBy(j) - gBy(j+1)) * ...
                        abs(gBz(k) - gBz(k+1));

                    % if this is the first timestep, find the particle that
                    % is closest to the center and log it's id
                    if startfilecheck == 1



                        % make a vector that is the central location,
                        % subtract from this a vector that is the x, y ,z
                        % particle positions, take magnitude, find smallest,
                        % that's the central particle
                        gridcenter = [ (gBy(j) + spacing/2 ), ...
                                 (gBz(k) + spacing/2)];

                        gridcenter_diffs = [y1,z1] - gridcenter;

                        min_gc_diff_arr = abs(gridcenter_diffs(:,1).*gridcenter_diffs(:,1) + ...
                            gridcenter_diffs(:,2).*gridcenter_diffs(:,2));

                        % find the minimum distance, this is the central
                        % particle which we'll consider the sensor
                        % mgdi = min_gc_diff_ind
                        [~,mgdi] = min(min_gc_diff_arr);

                        if j == s_i
                            store_psensor_ID(k,1) = p_id(mgdi); % this one is just logged initially
                            store_psensor_Fx{k,1} = [store_psensor_Fx{k,1}, fx(mgdi)];
                            store_psensor_Fy{k,1} = [store_psensor_Fy{k,1}, fy(mgdi)];
                            store_psensor_Fz{k,1} = [store_psensor_Fz{k,1}, fz(mgdi)];
                            store_psensor_x{k,1} = [store_psensor_x{k,1}, x1(mgdi)];
                            store_psensor_y{k,1} = [store_psensor_y{k,1}, y1(mgdi)];
                            store_psensor_z{k,1} = [store_psensor_z{k,1}, z1(mgdi)];
                            store_psensor_VY{k,1} = [store_psensor_VY{k,1}, vy(mgdi)];
                            store_psensor_VX{k,1} = [store_psensor_VX{k,1}, vx(mgdi)];
                            store_psensor_VZ{k,1} = [store_psensor_VZ{k,1}, vz(mgdi)];

                            store_psensor_stress_xx{k,1} = [store_psensor_stress_xx{k,1}, stress_xx(mgdi)];
                            store_psensor_stress_yy{k,1} = [store_psensor_stress_yy{k,1}, stress_yy(mgdi)];
                            store_psensor_stress_zz{k,1} = [store_psensor_stress_zz{k,1}, stress_zz(mgdi)];
                            store_psensor_stress_xy{k,1} = [store_psensor_stress_xy{k,1}, stress_xy(mgdi)];
                            store_psensor_stress_xz{k,1} = [store_psensor_stress_xz{k,1}, stress_xz(mgdi)];
                            store_psensor_stress_yz{k,1} = [store_psensor_stress_yz{k,1}, stress_yz(mgdi)];

                            store_psensor_delta{k,1} = [store_psensor_delta{k,1}, overlaps(p_id(mgdi))];
                        end




                    else % otherwise, just use the passed in particle ID
                        % look up particle ID

                        if j == s_i
                            this_particle_ind = find(p_id == store_psensor_ID(k,1));


                            if isempty(this_particle_ind)
                                % for this version  we don't care about the
                                % particle sensors, so just set arbitrarily
                                mgdi = 1;
                            else
                                mgdi = this_particle_ind;
                            end


                            store_psensor_Fx{k,1} = [store_psensor_Fx{k,1}, fx(mgdi)];
                            store_psensor_Fy{k,1} = [store_psensor_Fy{k,1}, fy(mgdi)];
                            store_psensor_Fz{k,1} = [store_psensor_Fz{k,1}, fz(mgdi)];
                            store_psensor_x{k,1} = [store_psensor_x{k,1}, x1(mgdi)];
                            store_psensor_y{k,1} = [store_psensor_y{k,1}, y1(mgdi)];
                            store_psensor_z{k,1} = [store_psensor_z{k,1}, z1(mgdi)];
                            store_psensor_VY{k,1} = [store_psensor_VY{k,1}, vy(mgdi)];
                            store_psensor_VX{k,1} = [store_psensor_VX{k,1}, vx(mgdi)];
                            store_psensor_VZ{k,1} = [store_psensor_VZ{k,1}, vz(mgdi)];

                            store_psensor_stress_xx{k,1} = [store_psensor_stress_xx{k,1}, stress_xx(mgdi)];
                            store_psensor_stress_yy{k,1} = [store_psensor_stress_yy{k,1}, stress_yy(mgdi)];
                            store_psensor_stress_zz{k,1} = [store_psensor_stress_zz{k,1}, stress_zz(mgdi)];
                            store_psensor_stress_xy{k,1} = [store_psensor_stress_xy{k,1}, stress_xy(mgdi)];
                            store_psensor_stress_xz{k,1} = [store_psensor_stress_xz{k,1}, stress_xz(mgdi)];
                            store_psensor_stress_yz{k,1} = [store_psensor_stress_yz{k,1}, stress_yz(mgdi)];

                            store_psensor_delta{k,1} = [store_psensor_delta{k,1}, overlaps(store_psensor_ID(k,1))];
                        end
                    end

    end
end

% bin the state data
[volumebins, Fmagbins, Fxbins, Fybins, Fzbins, VXbins, VYbins, VZbins, Vmagbins, ...
    delta_bins, bins_stress_xx, bins_stress_yy, bins_stress_zz, ...
    bins_stress_xy, bins_stress_xz, bins_stress_yz, bins_delta_total] ...
    = chunkParticles_overlaps(gBy, gBz, spacing, y1, z1, r1, vx, vy, vz, fx, fy, fz, m, evallim, ...
    overlaps, overlaps_total, stress_xx, stress_yy, stress_zz, stress_xy, stress_xz, stress_yz,p_id);

% get volume fraction of each cell
volfracbins = volumebins./gridvolbins;

vmag = (vx.*vx + vy.*vy + vz.*vz).^(1/2);

% voronoi tesselation
[vor_x, vor_y, vor_z, voro_vol,vor_r] = voronize(p_id,x1,y1,z1,r1,wall,zmax,x_lo,x_hi);


% bin the voronoi data
[voro_bins] = ...
    chunkParticles_voro(gBy, gBz, vor_y, vor_z, str2double(vor_r), voro_vol);


% find bed height (top sensor particles average height)
if startfilecheck == 1
    % get height of xy grid overlay and the particle IDS (these are the
    % "height sensors"
    [heights,bins_particle_ids] = findHeights_overlaps_init(z1,y1,gBy,r1,spacing,vmag,p_id,SF);

else
    [heights,bins_particle_ids] = findHeights_overlaps_allother(z1,y1,gBy,r1,spacing,p_id,bin_particles_tocheck,SF);

end


% if it is the last time step, use the bed heights to define a new
% voronization bed height. use max height in measurement region
voro_bins_end = 0;
if end_tstep_ind

    l_h = length(heights)
    h_start = round(SF*0.5/spacing) % find index of half the channl
    h_end = l_h - h_start

    % find height to use - height of particles more than 50 cm from a wall
    heights(h_start:h_end)
    height_end = mean(r1) + nanmean(heights(h_start:h_end))

    % define new grid bounds based on height
    % limits of evaluation region: [xlo, ylo, zlo; xhi, yhi, zhi]
    evallim_end = [x_lo,ystart,0;
                x_hi,(wall),height_end];

    % determine the grid y and z grid bounds. We average over x so there is only
    % 1 x index
    gBx_end = [evallim_end(1,1), evallim_end(2,1)];
    gBy_end = (evallim_end(1,2):spacing:evallim_end(2,2));
    gBz_end = (evallim_end(1,3):spacing:evallim_end(2,3))

    % reset the end value, this will stretch the last grid a bit, but will make
    % sure no area is left off
    gBy_end(end) = evallim_end(2,2);
    gBz_end(end) = evallim_end(2,3);


    numgridy_end = length(gBy_end) - 1;
    numgridz_end = length(gBz_end) - 1;

    gridvolbins_end = zeros(numgridy_end,numgridz_end);


    % perform the voronoi tesselation
    [~, vor_y_end, vor_z_end, voro_vol_end,vor_r_end] = voronize(p_id,x1,y1,z1,r1,wall,zmax,x_lo,x_hi);

    % perform the binning of the voronoi volumes
    [voro_bins_end] = ...
        chunkParticles_voro(gBy_end, gBz_end, vor_y_end, vor_z_end, str2double(vor_r_end), voro_vol_end);

end



end % for all y grids



% end function
