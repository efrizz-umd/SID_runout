% chunkParticles_voro.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [voro_bins] = ...
    chunkParticles_voro(gBy, gBz, y1, z1, r1, voro_vol)

% ************************************************************************
% This function performs the binning of particle voronoi volume (local packing fraction)
%
% inputs
% gBy, gBz - grid points in the y and z direction, respectively
% y1, z1 - particle positions (m)
% r1 - particle radius (m)
% voro_vol - particle voronoi volumes
%
% outputs
% - voro_bins - is the virtual sensor storage array for voronoi volume
% ************************************************************************

format long

numgridy = length(gBy) - 1;
numgridz = length(gBz) - 1;

numparticles = length(y1);

voro_bins = zeros(numgridy,numgridz);
bins_numparticles = zeros(numgridy,numgridz);

if numparticles > length(voro_vol)
    keyboard
end


% step through each particle from output. the particle looks up its
% position in the grid and places all into the grid bins based on orientation with grid bounds
% given our averaging approach, we ignore particles that intersect grid bounds to simplify the procedure
for k = 1:numparticles

    pcentery = y1(k);
    pcenterz = z1(k);
    r = r1(k);
    % voronoi volume
    vv = voro_vol(k);

    % find the index containing particle center
    diffy = (pcentery-gBy).';
    diffz = (pcenterz-gBz).';
    [~,mindexgby] = min(abs(diffy));
    [~,mindexgbz] = min(abs(diffz));
    % establish the particle gridbounds - closest boundaries to particle center
    pgby = gBy(mindexgby);
    pgbz = gBz(mindexgbz);
    % need to find orientation of particle in relation to gridbound (which way
    % to cut sphere)
    pdist(1,1) = pgby - pcentery;
    pdist(1,2) = pgbz - pcenterz;

    % first, determine which bin the mindex actually belogs to
    % this is not necessaryliy the same as the index returned for the closest
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
    %keyboard
    mindex(1,1) = length(find(pointerarry(:,1)));
    mindex(1,2) = length(find(pointerarrz(:,1)));

    % grab the particle overlap based on ID pointer - all other quantities
    % are output by particle so can just index those with k


    % check both distances. if greater than r, particle is contained wholly
    % in found grid
    if abs(pdist(1)) >= r && abs(pdist(2)) >= r
        if (mindex(1) <= numgridy) && (mindex(1) > 0) && (mindex(2) <= numgridz) && (mindex(2) > 0)
            voro_bins(mindex(1),mindex(2)) = voro_bins(mindex(1),mindex(2)) + volsphere(r)/vv;
            bins_numparticles(mindex(1),mindex(2)) = bins_numparticles(mindex(1),mindex(2)) + 1;

            %keyboard
        end
    end


end

% divide all by number of particles in the bin to get averages
voro_bins = voro_bins./bins_numparticles;



% keyboard
end
