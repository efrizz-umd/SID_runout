% findWindowAverages.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [F_rowavg,F_colavg,vel_rowavg,vel_colavg,vfpdiff_rowavg, vfpdiff_colavg, ...
F_rowstd,F_colstd,vel_rowstd,vel_colstd,vfpdiff_rowstd,vfpdiff_colstd, ...
F_rowavg_partial,F_colavg_partial,vel_rowavg_partial,vel_colavg_partial,vfpdiff_rowavg_partial,vfpdiff_colavg_partial, ...
F_rowstd_partial,F_colstd_partial,vel_rowstd_partial,vel_colstd_partial,vfpdiff_rowstd_partial,vfpdiff_colstd_partial, ...
F_rowavg_compact,F_colavg_compact,vel_rowavg_compact,vel_colavg_compact,vfpdiff_rowavg_compact,vfpdiff_colavg_compact, ...
F_rowstd_compact,F_colstd_compact,vel_rowstd_compact,vel_colstd_compact,vfpdiff_rowstd_compact,vfpdiff_colstd_compact] = ...
    findWindowAverages(ynumgrids,znumgrids,Ltp,store_vf_pdiff,store_vel,store_F, ssdilationind_full)

% ************************************************************************
% This function takes a storage array and computes the averages across columns and row
% for varying types (full dilation, partial, compaction). Inputs and outputs are similar
% to the other averaging functions, this one is not currently used.
%
% ************************************************************************


% *********** full dilation *********** %
%
%means
F_rowavg = zeros(znumgrids,Ltp);
F_colavg = zeros(ynumgrids,Ltp);
vel_rowavg = zeros(znumgrids,Ltp);
vel_colavg = zeros(ynumgrids,Ltp);
vfpdiff_rowavg = zeros(znumgrids,Ltp);
vfpdiff_colavg = zeros(ynumgrids,Ltp);
% standard deviations
F_rowstd = zeros(znumgrids,Ltp);
F_colstd = zeros(ynumgrids,Ltp);
vel_rowstd = zeros(znumgrids,Ltp);
vel_colstd = zeros(ynumgrids,Ltp);
vfpdiff_rowstd = zeros(znumgrids,Ltp);
vfpdiff_colstd = zeros(ynumgrids,Ltp);

% *********** partial dilation *********** %
%
%means
F_rowavg_partial = zeros(znumgrids,Ltp);
F_colavg_partial = zeros(ynumgrids,Ltp);
vel_rowavg_partial = zeros(znumgrids,Ltp);
vel_colavg_partial = zeros(ynumgrids,Ltp);
vfpdiff_rowavg_partial = zeros(znumgrids,Ltp);
vfpdiff_colavg_partial = zeros(ynumgrids,Ltp);
% standard deviations
F_rowstd_partial = zeros(znumgrids,Ltp);
F_colstd_partial = zeros(ynumgrids,Ltp);
vel_rowstd_partial = zeros(znumgrids,Ltp);
vel_colstd_partial = zeros(ynumgrids,Ltp);
vfpdiff_rowstd_partial = zeros(znumgrids,Ltp);
vfpdiff_colstd_partial = zeros(ynumgrids,Ltp);


% *********** compaction *********** %
%
%means
F_rowavg_compact = zeros(znumgrids,Ltp);
F_colavg_compact = zeros(ynumgrids,Ltp);
vel_rowavg_compact = zeros(znumgrids,Ltp);
vel_colavg_compact = zeros(ynumgrids,Ltp);
vfpdiff_rowavg_compact = zeros(znumgrids,Ltp);
vfpdiff_colavg_compact = zeros(ynumgrids,Ltp);
% standard deviations
F_rowstd_compact = zeros(znumgrids,Ltp);
F_colstd_compact = zeros(ynumgrids,Ltp);
vel_rowstd_compact = zeros(znumgrids,Ltp);
vel_colstd_compact = zeros(ynumgrids,Ltp);
vfpdiff_rowstd_compact = zeros(znumgrids,Ltp);
vfpdiff_colstd_compact = zeros(ynumgrids,Ltp);


% column average - based on dilation windows
for j = 1:ynumgrids
        % reset counter
        ccount = 1;
        ccount_partial = 1;
        ccount_compact = 1;
        for k = 1:znumgrids

            if ssdilationind_full(k,j) == 50

                if max(store_vf_pdiff{k,j}) < 20

                    tempc_vfpdiff(ccount,:) = store_vf_pdiff{k,j};
                    tempc_vel(ccount,:)  = store_vel{k,j};
                    tempc_F(ccount,:) = store_F{k,j};

                    % increment counter
                    ccount = ccount + 1;

                end
            elseif ssdilationind_full(k,j) == -50
                    tempc_vfpdiff_partial(ccount_partial,:) = store_vf_pdiff{k,j};
                    tempc_vel_partial(ccount_partial,:)  = store_vel{k,j};
                    tempc_F_partial(ccount_partial,:) = store_F{k,j};

                    % increment counter
                    ccount_partial = ccount_partial + 1;
            elseif ssdilationind_full(k,j) == -100
                    tempc_vfpdiff_compact(ccount_compact,:) = store_vf_pdiff{k,j};
                    tempc_vel_compact(ccount_compact,:)  = store_vel{k,j};
                    tempc_F_compact(ccount_compact,:) = store_F{k,j};

                    % increment counter
                    ccount_compact = ccount_compact + 1;

            end



        end

        % full dilation
        if ccount ~= 1

            if ccount == 2
                vfpdiff_colavg(j,:) = tempc_vfpdiff;
                vel_colavg(j,:) = tempc_vel;
                F_colavg(j,:) = tempc_F;

                vfpdiff_colstd(j,:) = 0;
                vel_colstd(j,:) = 0;
                F_colstd(j,:) = 0;
            else
                vfpdiff_colavg(j,:) = nanmean(tempc_vfpdiff);
                vel_colavg(j,:) = nanmean(tempc_vel);
                F_colavg(j,:) = nanmean(tempc_F);

                vfpdiff_colstd(j,:) = nanstd(tempc_vfpdiff);
                vel_colstd(j,:) = nanstd(tempc_vel);
                F_colstd(j,:) = nanstd(tempc_F);
            end
            % clear temporary variables
            clear tempc_vfpdiff tempc_vel tempc_F
        end

        % partial dilation
        if ccount_partial ~= 1

            if ccount_partial == 2
                vfpdiff_colavg_partial(j,:) = tempc_vfpdiff_partial;
                vel_colavg_partial(j,:) = tempc_vel_partial;
                F_colavg_partial(j,:) = tempc_F_partial;

                vfpdiff_colstd_partial(j,:) = 0;
                vel_colstd_partial(j,:) = 0;
                F_colstd_partial(j,:) = 0;
            else
                vfpdiff_colavg_partial(j,:) = nanmean(tempc_vfpdiff_partial);
                vel_colavg_partial(j,:) = nanmean(tempc_vel_partial);
                F_colavg_partial(j,:) = nanmean(tempc_F_partial);

                vfpdiff_colstd_partial(j,:) = nanstd(tempc_vfpdiff_partial);
                vel_colstd_partial(j,:) = nanstd(tempc_vel_partial);
                F_colstd_partial(j,:) = nanstd(tempc_F_partial);
            end
            % clear temporary variables
            clear tempc_vfpdiff_partial tempc_vel_partial tempc_F_partial
        end

        % compaction
        if ccount_compact ~= 1

            if ccount_compact == 2

                vfpdiff_colavg_compact(j,:) = tempc_vfpdiff_compact;
                vel_colavg_compact(j,:) = tempc_vel_compact;
                F_colavg_compact(j,:) = tempc_F_compact;

                vfpdiff_colstd_compact(j,:) = 0;
                vel_colstd_compact(j,:) = 0;
                F_colstd_compact(j,:) = 0;
            else
                vfpdiff_colavg_compact(j,:) = nanmean(tempc_vfpdiff_compact);
                vel_colavg_compact(j,:) = nanmean(tempc_vel_compact);
                F_colavg_compact(j,:) = nanmean(tempc_F_compact);

                vfpdiff_colstd_compact(j,:) = nanstd(tempc_vfpdiff_compact);
                vel_colstd_compact(j,:) = nanstd(tempc_vel_compact);
                F_colstd_compact(j,:) = nanstd(tempc_F_compact);
            end
            % clear temporary variables
            clear tempc_vfpdiff_compact tempc_vel_compact tempc_F_compact
        end

end

% row average - based on dilation windows
for k = 1:znumgrids
        % reset counter
        rcount = 1;
        rcount_partial = 1;
        rcount_compact = 1;
        for j = 1:ynumgrids

            % dilation
            if ssdilationind_full(k,j) == 50

                if max(store_vf_pdiff{k,j}) < 20
                    tempr_vfpdiff(rcount,:) = store_vf_pdiff{k,j};
                    tempr_vel(rcount,:)  = store_vel{k,j};
                    tempr_F(rcount,:) = store_F{k,j};

                    % increment counter
                    rcount = rcount + 1;
                end
            end

            % partial dilation
            if ssdilationind_full(k,j) == -50

                tempr_vfpdiff_partial(rcount_partial,:) = store_vf_pdiff{k,j};
                tempr_vel_partial(rcount_partial,:)  = store_vel{k,j};
                tempr_F_partial(rcount_partial,:) = store_F{k,j};

                % increment counter
                rcount_partial = rcount_partial + 1;

            end

            % compaction
            if ssdilationind_full(k,j) == -100

                tempr_vfpdiff_compact(rcount_compact,:) = store_vf_pdiff{k,j};
                tempr_vel_compact(rcount_compact,:)  = store_vel{k,j};
                tempr_F_compact(rcount_compact,:) = store_F{k,j};

                % increment counter
                rcount_compact = rcount_compact + 1;

            end

        end


        % dilation
        if rcount ~= 1
            % average
            if rcount == 2
                vfpdiff_rowavg(k,:) = tempr_vfpdiff;
                vel_rowavg(k,:) = tempr_vel;
                F_rowavg(k,:) = tempr_F;

                vfpdiff_rowstd(k,:) = 0;
                vel_rowstd(k,:) = 0;
                F_rowstd(k,:) = 0;
            else
                vfpdiff_rowavg(k,:) = nanmean(tempr_vfpdiff);
                vel_rowavg(k,:) = nanmean(tempr_vel);
                F_rowavg(k,:) = nanmean(tempr_F);

                vfpdiff_rowstd(k,:) = nanstd(tempr_vfpdiff);
                vel_rowstd(k,:) = nanstd(tempr_vel);
                F_rowstd(k,:) = nanstd(tempr_F);
            end
            % clear temporary variables
            clear tempr_vfpdiff tempr_vel tempr_F
        end

        % partial dilation
        if rcount_partial ~= 1
            % average
            if rcount_partial == 2
                vfpdiff_rowavg_partial(k,:) = tempr_vfpdiff_partial;
                vel_rowavg_partial(k,:) = tempr_vel_partial;
                F_rowavg_partial(k,:) = tempr_F_partial;

                vfpdiff_rowstd_partial(k,:) = 0;
                vel_rowstd_partial(k,:) = 0;
                F_rowstd_partial(k,:) = 0;
            else
                vfpdiff_rowavg_partial(k,:) = nanmean(tempr_vfpdiff_partial);
                vel_rowavg_partial(k,:) = nanmean(tempr_vel_partial);
                F_rowavg_partial(k,:) = nanmean(tempr_F_partial);

                vfpdiff_rowstd_partial(k,:) = nanstd(tempr_vfpdiff_partial);
                vel_rowstd_partial(k,:) = nanstd(tempr_vel_partial);
                F_rowstd_partial(k,:) = nanstd(tempr_F_partial);
            end
            % clear temporary variables
            clear tempr_vfpdiff_partial tempr_vel_partial tempr_F_partial
        end


        % compaction
        if rcount_compact ~= 1
            % average
            if rcount_compact == 2
                vfpdiff_rowavg_compact(k,:) = tempr_vfpdiff_compact;
                vel_rowavg_compact(k,:) = tempr_vel_compact;
                F_rowavg_compact(k,:) = tempr_F_compact;

                vfpdiff_rowstd_compact(k,:) = 0;
                vel_rowstd_compact(k,:) = 0;
                F_rowstd_compact(k,:) = 0;
            else
                vfpdiff_rowavg_compact(k,:) = nanmean(tempr_vfpdiff_compact);
                vel_rowavg_compact(k,:) = nanmean(tempr_vel_compact);
                F_rowavg_compact(k,:) = nanmean(tempr_F_compact);

                vfpdiff_rowstd_compact(k,:) = nanstd(tempr_vfpdiff_compact);
                vel_rowstd_compact(k,:) = nanstd(tempr_vel_compact);
                F_rowstd_compact(k,:) = nanstd(tempr_F_compact);
            end
            % clear temporary variables
            clear tempr_vfpdiff_compact tempr_vel_compact tempr_F_compact
        end
end


end
