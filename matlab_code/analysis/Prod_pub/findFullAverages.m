% findFullaverages.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [quantity_rowavg_full,quantity_colavg_full,quantity_rowstd_full,quantity_colstd_full] = ...
    findFullAverages(ynumgrids,znumgrids,Ltp,store)

% ************************************************************************
% This function takes a storage array and computes the averages across columns and row,
% for a single time step. Operates on cell arrays
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
% - Ltp - numer of time steps
% ************************************************************************

%means
quantity_rowavg_full = zeros(znumgrids,Ltp);
quantity_colavg_full = zeros(ynumgrids,Ltp);

% standard deviations
quantity_rowstd_full = zeros(znumgrids,Ltp);
quantity_colstd_full = zeros(ynumgrids,Ltp);


% column average - based on all windows
for j = 1:ynumgrids
    for k = 1:znumgrids
        tempc_quantity_full(k,:) = store{k,j};
    end

    quantity_colavg_full(j,:) = nanmean(tempc_quantity_full);

    quantity_colstd_full(j,:) = nanstd(tempc_quantity_full);


    % clear temporary variables
    clear tempc_quantity_full

end

% row average - based on all windows
for k = 1:znumgrids
    for j = 1:ynumgrids
        tempr_quantity_full(j,:) = store{k,j};
    end

    quantity_rowavg_full(k,:) = nanmean(tempr_quantity_full);

    quantity_rowstd_full(k,:) = nanstd(tempr_quantity_full);


    % clear temporary variables
    clear tempr_quantity_full

end


end
