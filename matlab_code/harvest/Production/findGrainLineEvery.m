function [grainline, channellength] ...
               = findGrainLineEvery(filename,ii,SF)

% ************************************************************************
% This function takes a LIGGGHTS post file and plots the density of the
% granular material inside the simulation domain as a function of distance
% from the wall. The density is calculated by dividing the domain into a
% number of different bins of a given "spacing"(input 3)
%
% input
% file - file location
%
% output
% grainline is the average max height of the particle
% ************************************************************************

%temp = strcat('dump',strcat(num2str(step),'.post'));
%tempo = strcat('dumpcad',strcat(num2str(step),'.stl'));
%filename = strcat(file,temp);
%stlname = strcat(file,tempo);

a = strcat('dump',strcat(num2str(ii),'.post'));
filezero = strcat(filename,a);

 
%% Find bounds
% delimiter = ' ';
% startRow = 6;
% endRow = 8;
% formatSpec = '%f%f%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';
% fileID = fopen(filename,'r');
% %keyboard
% % import simulation bounds from the header
% dataArray = textscan(fileID, formatSpec, endRow-startRow+1, ...
% 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow-1, 'ReturnOnError', false);
% fclose(fileID);
% 
% % Allocate imported array to column variable names. hiBound is the side
% % with the wall, and hiBound(2) is the wall-dependent boundary
% loBound = dataArray{:, 1};
% hiBound = dataArray{:, 2};
% 
% 
% clearvars delimiter startRow formatSpec fileID dataArray ans;
%% Read Positions

delimiter = ' ';
startRow = 10;

% asterisk after the delimiter marker (%) means skip that field
%formatSpec = '%*s%f%f%f%f%*s%*s%*s%f%*s%*s%[^\n\r]';
formatSpec = '%*s%*s%f%f%f%*s%*s%*s%f%*s%*s%*s%*s%[^\n\r]';

fileID = fopen(filezero,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);
%id = dataArray{:, 1};
x1 = dataArray{:, 1};
y1 = dataArray{:, 2};
z1 = dataArray{:, 3};
%vx = dataArray{:, 5};
% vy = dataArray{:, 6};
% vz = dataArray{:, 7};
 r1 = dataArray{:, 4};
% fx = dataArray{:, 9};
% fy = dataArray{:, 10};
% fz = dataArray{:, 11};
% m = dataArray{:, 12};

numparticles = length(x1);

% determine total surface area
xmin = min(x1);
xmax = max(x1);
ymax = max(y1);

%channellength = round(ymax);
channellength = ymax;

% setup a surface area grid
spacing = SF*.01; % 1 cm

gbx = 0:spacing:(SF*.02);
gby = 0:spacing:ymax;

numgridx = length(gbx) - 1;
numgridy = length(gby) - 1;

% data storage for maximum heights per bin
maxheight = zeros(numgridx,numgridy);



for k = 1:numparticles
    
    pcenterx = x1(k) + abs(xmin);
    pcentery = y1(k);
    r = r1(k);
    
    % find the index containing particle center
    diffx = (pcenterx-gbx).';
    diffy = (pcentery-gby).';
    
    
%     [mindisty,mindexgby] = min(abs(diffy));
%     [mindistz,mindexgbz] = min(abs(diffz));
%     % establish the particle gridbounds - closest boundaries to particle center
%     pgby = gBy(mindexgby); 
%     pgbz = gBz(mindexgbz); 
%     % need to find orientation of particle in relation to gridbound (which way
%     % to cut sphere)
%     pdist(1,1) = pgby - pcentery;
%     pdist(1,2) = pgbz - pcenterz;

    % first, determine which bin the mindex actually belogs to
    % this is not necessaryliy the same as the index returned for the closest
    % grid bound
    % make normal vectors  (1D) vectors which point from the plane defined by
    % the gridbound to +y or from 
    normalvec = [100]; 
    %duplicate diffs,s then overwrite
    diffdirx = diffx;
    diffdiry = diffy;
    for j = 1:length(diffx)
        diffdirx(j) = sign(dot(diffx(j),normalvec));
    end
    for j = 1:length(diffy)
        diffdiry(j) = sign(dot(diffy(j),normalvec));
    end
    % subtract one from the results, this makes the grid that contains the
    % particle bounded by 0,-2 instead of 1,-1, then take the ~ and find
    % the number first instance of 1
    diffdirx = diffdirx - 1;
    diffdiry = diffdiry - 1;
    pointerarrx = ~diffdirx;
    pointerarry = ~diffdiry;
    %keyboard
    mindex(1,1) = length(find(pointerarrx(:,1)));
    mindex(1,2) = length(find(pointerarry(:,1)));
    
    pheight = z1(k) + r;
%     if k == 1447429
%         keyboard
%     end
    
    
    if mindex(1,1) <= numgridx && mindex(1,2) <= numgridy && mindex(1,1) > 0 && mindex(1,2) > 0
        if pheight > maxheight(mindex(1,1),mindex(1,2))
            maxheight(mindex(1,1),mindex(1,2)) = pheight;
        end
    end
    %keyboard
end

grainline = mean(maxheight(:));



% radius = mean(r1);
% xlen = max(x1) - min(x1);
% ylen = max(y1);
% Asurf = xlen*ylen;
% numkeep = floor(Asurf/(pi*radius^2));
% grainline = mean(maxk(z1,numkeep));



