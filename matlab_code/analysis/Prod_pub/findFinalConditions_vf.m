% findFinalConditions_vf.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [rowavg_init, rowstd_init, colavg_init, colstd_init] = ...
    findFinalConditions_vf(store,ynumgrids,znumgrids)

% ************************************************************************
% This function takes a storage array and computes the averages across columns and row,
% for a single time step
%
% % ----------- output ----------- %
% - rowavg_init - average across the rows (average at a given depth)
% - rowstd_init - std across the rows (std. dev. at a given depth)
% - colavg_init - average across the columns (average at given radial position)
% - colstd_init - std across the columns (std. dev. at a given radial position)

% % ----------- intput ----------- %
% - store - a double array of some values you want to average
% - ynumgrids - number of radial grids
% - znumgrids - number of depth grids
% ************************************************************************

%% initialize storage
%means
rowavg_init = NaN(1,znumgrids);
colavg_init = NaN(1,ynumgrids);

% standard deviations
rowstd_init = NaN(1,znumgrids);
colstd_init = NaN(1,ynumgrids);


%% loop and find by row or by colum

% the rows
for j = 1:znumgrids
    tempR = NaN(1,ynumgrids);

    % build the initial conditions across the row
    for k = 1:ynumgrids
        tempR(1,k) = store(j,k);

        if store(j,k) <= 0.001
            tempR(1,k) = NaN;
        end
    end

    % mean and standard deviations
    rowavg_init(1,j) = nanmean(tempR);
    rowstd_init(1,j) = nanstd(tempR);


end

% the columns
for j = 1:ynumgrids
    tempC = NaN(1,znumgrids);

    % build the initial conditions across the columns
    for k = 1:znumgrids
        tempC(1,k) = store(k,j);

        if store(k,j) <= 0.001
            tempC(1,k) = NaN;
        end
    end

    % mean and standard deviations
    colavg_init(1,j) = nanmean(tempC);
    colstd_init(1,j) = nanstd(tempC);

end


% zero protection
%means
rowavg_init(rowavg_init == 0) = NaN;
colavg_init(colavg_init == 0) = NaN;

% standard deviations
rowstd_init(rowstd_init == 0) = NaN;
colstd_init(colstd_init == 0) = NaN;


end
