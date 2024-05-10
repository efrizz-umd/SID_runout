% chunkParticles_overlaps.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [volumebins, bins_totalF, bins_Fx, bins_Fy, bins_Fz, ...
    bins_vx, bins_vy, bins_vz, bins_totalV,bins_delta, ...
    bins_stress_xx, bins_stress_yy, bins_stress_zz, ...
    bins_stress_xy, bins_stress_xz, bins_stress_yz, bins_delta_total] = ...
    chunkParticles_overlaps(gBy, gBz, spacing, y1, z1, r1, vx, vy, vz, fx, fy, fz, m, evallim, ...
    delta,delta_total, stress_xx, stress_yy, stress_zz, stress_xy, stress_xz, stress_yz,p_id)

% ************************************************************************
% This function performs the binning of particle state information (velocity,
% force, stress, overlap)
%
% inputs
% gBy, gBz - grid points in the y and z direction, respectively
% x1, y1, z1 - particle positions (m)
% r1 - particle radius (m)
% vx, vy, vz - particle velocities (m/s)
% fx, fy, fz - particle forces (N)
% evallim - if there are any filters or limiters to apply
% delta, delta_total - arrays of the average and total overlap per particle
% stress_** - stress in a given direction (LIGGGHTS stress units)
% p_id - particle IDs
%
% outputs
% - bins_* - is the virtual sensor storage array for a given quantity (F = force
% V = velocity, stress = particle stresses, voro = bins of packing fraction computed
% via voronoi tesselation). There will either be an associated direction (x, y, z)
% or magnitude ('mag') flag for each named particle.
% ************************************************************************

format long

% form array
%gB = [ymin:sy:ymax; zmin:sz:zmax];

% particle [y z]'
%pcenter = [9.99; 9.99];
%r = 2;

numgridy = length(gBy) - 1;
numgridz = length(gBz) - 1;

numparticles = length(y1);

% initialize storage
volumebins = zeros(numgridy,numgridz);
bins_numparticles = zeros(numgridy,numgridz);
bins_numparticles_full = zeros(numgridy,numgridz);
bins_totalF  = zeros(numgridy,numgridz);
bins_Fx = zeros(numgridy,numgridz);
bins_Fy = zeros(numgridy,numgridz);
bins_Fz = zeros(numgridy,numgridz);
bins_vx = zeros(numgridy,numgridz);
bins_vy = zeros(numgridy,numgridz);
bins_vz = zeros(numgridy,numgridz);
bins_totalV = zeros(numgridy,numgridz);
bins_mass = zeros(numgridy,numgridz);
bins_delta = zeros(numgridy,numgridz);
bins_delta_total = zeros(numgridy,numgridz);
bins_stress_xx = zeros(numgridy,numgridz);
bins_stress_yy = zeros(numgridy,numgridz);
bins_stress_zz = zeros(numgridy,numgridz);
bins_stress_xy = zeros(numgridy,numgridz);
bins_stress_xz = zeros(numgridy,numgridz);
bins_stress_yz = zeros(numgridy,numgridz);

% step through each particle from output. the particle looks up its
% position in the grid and places all into the grid bins based on orientation with grid bounds
% given our averaging approach, we ignore particles that intersect grid bounds to simplify the procedure
for k = 1:numparticles

    pcentery = y1(k);
    pcenterz = z1(k);
    r = r1(k);

    % find the index containing particle center
    diffy = (pcentery-gBy).';
    diffz = (pcenterz-gBz).';
    [mindisty,mindexgby] = min(abs(diffy));
    [mindistz,mindexgbz] = min(abs(diffz));
    % establish the particle gridbounds - closest boundaries to particle center
    pgby = gBy(mindexgby);
    pgbz = gBz(mindexgbz);
    % need to find orientation of particle in relation to gridbound
    pdist(1,1) = pgby - pcentery;
    pdist(1,2) = pgbz - pcenterz;

    % first, determine which bin the mindex actually belogs to
    % this is not necessarily the same as the index returned for the closest
    % grid bound
    % make normal vectors  (1D) vectors which point from the plane defined by
    % the gridbound to +y or from
    normalvec = [100];
    %duplicate diffs,s then overwrite
    diffdiry = diffy;
    diffdirz = diffz;
    for j = 1:length(diffy)
        diffdiry(j) = sign(dot(diffy(j),normalvec));
    end
    for j = 1:length(diffz)
        diffdirz(j) = sign(dot(diffz(j),normalvec));
    end
    % subtract one from the results, this makes the grid that contains the
    % particle bounded by 0,-2 instead of 1,-1, then take the ~ and find
    % the number first instance of 1
    diffdiry = diffdiry - 1;
    diffdirz = diffdirz - 1;
    pointerarry = ~diffdiry;
    pointerarrz = ~diffdirz;

    mindex(1,1) = length(find(pointerarry(:,1)));
    mindex(1,2) = length(find(pointerarrz(:,1)));

    % grab the particle overlap based on ID pointer - all other quantities
    % are output by particle so can just index those with k
    p_overlap = delta(p_id(k));
    p_overlap_total = delta_total(p_id(k));
    if isnan(p_overlap)
        p_overlap = 0;
    end
    if isnan(p_overlap_total)
        p_overlap_total = 0;
    end

    % first, check both distances. if greater than r, particle is contained wholly
    % in found grid
    if abs(pdist(1)) >= r && abs(pdist(2)) >= r
        if (mindex(1) <= numgridy) && (mindex(1) > 0) && (mindex(2) <= numgridz) && (mindex(2) > 0)
            volumebins(mindex(1),mindex(2)) = volumebins(mindex(1),mindex(2)) + volsphere(r);
            bins_numparticles(mindex(1),mindex(2)) = bins_numparticles(mindex(1),mindex(2)) + 1;
            bins_totalF(mindex(1),mindex(2))  = bins_totalF(mindex(1),mindex(2)) + norm([fx(k) fy(k) fz(k)]);
            bins_Fx(mindex(1),mindex(2)) = bins_Fx(mindex(1),mindex(2)) + (fx(k));
            bins_Fy(mindex(1),mindex(2)) = bins_Fy(mindex(1),mindex(2)) + (fy(k));
            bins_Fz(mindex(1),mindex(2)) = bins_Fz(mindex(1),mindex(2)) + (fz(k));
            bins_vx(mindex(1),mindex(2)) = bins_vx(mindex(1),mindex(2)) + (vx(k));
            bins_vy(mindex(1),mindex(2)) = bins_vy(mindex(1),mindex(2)) + (vy(k));
            bins_vz(mindex(1),mindex(2)) = bins_vz(mindex(1),mindex(2)) + (vz(k));
            bins_totalV(mindex(1),mindex(2)) = bins_totalV(mindex(1),mindex(2)) + norm([vx(k) vy(k) vz(k)]);
            bins_mass(mindex(1),mindex(2)) = bins_mass(mindex(1),mindex(2)) + m(k);

            bins_delta(mindex(1),mindex(2)) = bins_delta(mindex(1),mindex(2)) + p_overlap;
            bins_delta_total(mindex(1),mindex(2)) = bins_delta_total(mindex(1),mindex(2)) + p_overlap_total;

            bins_stress_xx(mindex(1),mindex(2)) = bins_stress_xx(mindex(1),mindex(2)) + stress_xx(k);
            bins_stress_yy(mindex(1),mindex(2)) = bins_stress_yy(mindex(1),mindex(2)) + stress_yy(k);
            bins_stress_zz(mindex(1),mindex(2)) = bins_stress_zz(mindex(1),mindex(2)) + stress_zz(k);
            bins_stress_xy(mindex(1),mindex(2)) = bins_stress_xy(mindex(1),mindex(2)) + stress_xy(k);
            bins_stress_xz(mindex(1),mindex(2)) = bins_stress_xz(mindex(1),mindex(2)) + stress_xz(k);
            bins_stress_yz(mindex(1),mindex(2)) = bins_stress_yz(mindex(1),mindex(2)) + stress_yz(k);

            bins_numparticles_full(mindex(1),mindex(2)) = bins_numparticles_full(mindex(1),mindex(2)) + 1;

        end
    end

end

% divide all by number of particles in the bin to get averages
bins_totalF  = bins_totalF./bins_numparticles_full;
bins_Fx = bins_Fx./bins_numparticles_full;
bins_Fy = bins_Fy./bins_numparticles_full;
bins_Fz = bins_Fz./bins_numparticles_full;
bins_vx = bins_vx./bins_numparticles_full;
bins_vy = bins_vy./bins_numparticles_full;
bins_vz = bins_vz./bins_numparticles_full;
bins_totalV = bins_totalV./bins_numparticles;
bins_delta = bins_delta./bins_numparticles_full;
bins_stress_xx = bins_stress_xx./bins_numparticles_full;
bins_stress_yy = bins_stress_yy./bins_numparticles_full;
bins_stress_zz = bins_stress_zz./bins_numparticles_full;
bins_stress_xy = bins_stress_xy./bins_numparticles_full;
bins_stress_xz = bins_stress_xz./bins_numparticles_full;
bins_stress_yz = bins_stress_yz./bins_numparticles_full;


% keyboard
end
