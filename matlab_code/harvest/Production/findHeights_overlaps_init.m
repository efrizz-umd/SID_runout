function [bins_average_height,bins_particle_ids_top] = findHeights_overlaps_init(z1,y1,gBy,r1,spacing, ...
    startfilecheck,p_id,SF)
%FINDHEIGHTS Summary of this function goes here
%   Detailed explanation goes here

% if this is the initial use, log the particle IDs. if not, just find the
% height of the input particle ids

numgridy = length(gBy) - 1;
numparticles = length(y1);
bins_y = cell(1,numgridy);
bins_r = cell(1,numgridy);
bins_particle_ids = cell(1,numgridy);
bins_particle_ids_top = cell(1,numgridy);
bins_average_height = zeros(1,numgridy);

mingBy = gBy(1);
maxgBy = gBy(end);


for k = 1:numparticles
    
    pcentery = y1(k);
    pcenterz = z1(k);
    r = r1(k);
    p_id_particle = p_id(k);
%     
%     if p_id(k) == 141914
%         keyboard
%     end
    
    if pcentery > mingBy && pcentery < maxgBy
        % find the index containing particle center
        diffy = (pcentery-gBy).';
        [~,mindexgby] = min(abs(diffy));
        % establish the particle gridbounds - closest boundaries to particle center
        pgby = gBy(mindexgby); 
        % need to find orientation of particle in relation to gridbound (which way
        % to cut sphere)
        pdist(1,1) = pgby - pcentery;

        % first, determine which bin the mindex actually belogs to
        % this is not necessaryliy the same as the index returned for the closest
        % grid bound
        % make normal vectors  (1D) vectors which point from the plane defined by
        % the gridbound to +y or from 
        normalvec = [100]; 
        %duplicate diffs,s then overwrite
        diffdiry = diffy;
        for j = 1:length(diffy)
            diffdiry(j) = sign(dot(diffy(j),normalvec));
        end
        % subtract one from the results, this makes the grid that contains the
        % particle bounded by 0,-2 instead of 1,-1, then take the ~ and find
        % the number first instance of 1
        diffdiry = diffdiry - 1;
        pointerarry = ~diffdiry;
        %keyboard
        mindex = length(find(pointerarry(:,1)));


        % log heights of particles in this index
        
        bins_y{1,mindex} = [bins_y{1,mindex}, pcenterz+r];
        bins_r{1,mindex} = [bins_r{1,mindex}, r];
        bins_particle_ids{1,mindex} = [bins_particle_ids{1,mindex}, p_id_particle];
    end
end

    % number of particles that would make up the surface
    mean_rad = mean(r1);
    mean_surf = pi*mean_rad*mean_rad;
    grid_area = spacing*0.02*SF; % hardcoded x dimension of 20 mm
    num_max = floor(grid_area/mean_surf);
    
    velthresh = 0.2; % 0.2 m/s, heuristically chosen based on ovito velocity filter
    
for jj = 1:numgridy
    
    % set particles with greater than 0.2 m/s velocity to 0 so they don't
    % get included. these particles are still flying and we wont consider
    % them. will probably chop off some particles from the height near
    % impact and early in the sim, but we don't care about those ^\O.o/^
%     the_y_bin = bins_y{1,jj};
%     the_y_bin(vmag >= velthresh) = 0;
    
%     if jj == 121
%         keyboard
%     end
    % average height of the heighest #num_max particles
    [p_heights,pid_thisgrid] = maxk(bins_y{1,jj},num_max);
    bins_average_height(1,jj) = mean(p_heights);
%     bins_particle_ids_top{1,jj} = p_id(bins_particle_ids{1,jj}(pid_thisgrid));
    bins_particle_ids_top{1,jj} = bins_particle_ids{1,jj}(pid_thisgrid);
    %keyboard
end

