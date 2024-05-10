% findFullMax.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [max_q_radial_fullw, max_qtime_radial_fullw, max_qdist_radial_fullw, max_qstd_radial_fullw , ...
    max_q_depth_fullw, max_qtime_depth_fullw, max_qdepth_fullw] ...
     = findFullMax(ynumgrids,znumgrids,spmtind, wspacing, windowrange,timeplot, ...
     q_colavg_full,q_colstd_full,q_rowavg_full)

 % ************************************************************************
 % This function takes a storage array and computes the averages across columns and row,
 % for a single time step. Operates on cell arrays. <DEPRECATED>, no longer used

 %
 % % ----------- output ----------- %
 % - quantity_rowavg_full - average across the rows (average at a given depth)
 % - quantity_rowstd_full - std across the rows (std. dev. at a given depth)
 % - quantity_colavg_full - average across the columns (average at given radial position)
 % - quantity_colstd_full - std across the columns (std. dev. at a given radial position)

 % % ----------- intput ----------- %
 % - store - a cell array of some values you want to average
 % - ynumgrids - number of radial grids
 % - znumgrids - number of depth grids
 % - timeplot - time array
 % - spmtind - shock period max time index
 % ************************************************************************

% ******** radial (depth averaged) *********
max_q_radial_fullw = zeros(1, ynumgrids);
max_qtime_radial_fullw = zeros(1, ynumgrids);
max_qdist_radial_fullw = zeros(1, ynumgrids);
max_qstd_radial_fullw = zeros(1, ynumgrids);

% limit findings to shock period only (1:spmtind)
for j=1:ynumgrids

    chbeg = (j-1)*wspacing + windowrange(1);
    chend = (j)*wspacing + windowrange(1);

    % find max and log time

    if spmtind > length(q_colavg_full(j,:))
        spmtind = length(q_colavg_full(j,:));
    end
    [max_q_radial_fullw(1,j), maxdex] = max(q_colavg_full(j,1:spmtind));
    max_qtime_radial_fullw(1,j) = timeplot(maxdex);
    max_qdist_radial_fullw(1,j) = mean([chbeg chend]);
    max_qstd_radial_fullw(1,j) = q_colstd_full(j,maxdex);

end

% ******** depth *********
% force
max_q_depth_fullw = zeros(1, znumgrids);
max_qtime_depth_fullw = zeros(1, znumgrids);
max_qdepth_fullw = zeros(1, znumgrids);

for j=1:znumgrids

    depthbeg = (j-1)*wspacing;
    depthend = (j)*wspacing;

    % find max and log time
    [max_q_depth_fullw(1,j), maxdex] = max(q_rowavg_full(j,1:spmtind));
    max_qtime_depth_fullw(1,j) = timeplot(maxdex);
    max_qdepth_fullw(1,j) = mean([depthbeg depthend]);



end

end
