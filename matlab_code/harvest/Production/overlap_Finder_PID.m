% overlap_Finder_PID.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [average_delta, total_delta] ...
               = overlap_Finder_PID(filename,numparticles)

% ************************************************************************
% This function finds the average and total overlap from a LIGGGHTS compute file
% (ex: dump_computes0.post). The output contains the overlap data for each particle in the system,
% sorted by particle ID
%
% inputs
% filename - file name corresponding to overlaps (ex: dump_computes0.post)
% num particles - total number of particle in the sim
%
% outputs 
% average_delta - array of each particle's average overlap, sorted by ID
% total_delta - array of each particle's total overlap, sorted by ID
% ************************************************************************


% file name
fn = "/Users/hartzell_lab/DOCUMENTS/COLDSPOTS/MATLAB ANALYSIS/ShockInducedDilation/PlusOverlapAndStress/post_computes/dump_computes0.post";

delimiter = ' ';
startRow = 10;

% asterisk after the delimiter marker (%) means skip that field
formatSpec = '%f%f%f%[^\n\r]';

% open and read data
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
fclose(fileID);


% assign extracted data to arrays
id1 = dataArray{:, 1};
id2 = dataArray{:, 2};
deltas = dataArray{:, 3};

L_entries = length(id1);

ids = [id1; id2];
unique_ids = unique(ids);
num_entries = max(unique_ids);
L_ids = length(unique_ids);

% length of storage arrays will have the same length as number of
% particles, so that the particle ID can be used in indexing
store_overlaps = cell(1,numparticles);
store_numcontacts = cell(1,numparticles);

% final array for storing mean overlap
average_delta = zeros(1,numparticles);
total_delta = zeros(1,numparticles);

% accumulate overlaps per particle
for jj = 1:L_entries

    temp_id1 = id1(jj,1);
    temp_id2 = id2(jj,1);

    store_overlaps{1,temp_id1} = [store_overlaps{1,temp_id1} deltas(jj,1)];
    store_overlaps{1,temp_id2} = [store_overlaps{1,temp_id2} deltas(jj,1)];

end

% sum or average each particles overlap state
for jj = 1:num_entries
    average_delta(1,jj) = mean(store_overlaps{1,jj});
    total_delta(1,jj) = sum(store_overlaps{1,jj});
end
