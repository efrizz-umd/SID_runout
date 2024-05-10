% HarvestAndPlotFun_production_loadonly.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [grainstring, channelstring, porstring, ...
    channel, vy_max, vy_maxtime, percentofmaxvel, Fmag_max, percentofmaxF, ...
    Fmag_colavg_init, Fmag_colstd_init, vy_colavg_init, vy_colstd_init, vf_colavg_init, vf_colstd_init, ...
    depthplot,Fmag_rowavg_init,Fmag_rowstd_init, vy_rowavg_init, vy_rowstd_init, vf_rowavg_init, vf_rowstd_init, ...
    max_Fmagdist_radial_fullw,max_Fmag_radial_fullw,max_Fmagstd_radial_fullw, max_Fmagtime_radial_fullw, ...
    max_vmagdist_radial_fullw,max_vmag_radial_fullw,max_vmagstd_radial_fullw, max_vmagtime_radial_fullw, ...
    shockmax,shockmaxstd, shockmin,shockminstd,shockmintime,shock_runout_V, shock_speed_V, ...
    shock_runout_F,shock_speed_F,shockmin_F,shockspeedslope, ...dilationpdiff_SS_tp_f_a,dilationpdiff_SS_tp_p_a,...
    timeplot,regionstartdepth,regionstopdepth,Fmag_logref_mean,Fmag_logref_std, ...
    bulk_dilation_vsdepth,bulk_dilation_stddv_vsdepth,F_shocktrack_mean,shockspeedslope_arrival, ...
    bulk_dilation_vsradial, bulk_dilation_stddv_vsradial, shockspeedslope_maxforce, ...
    height_pdiff_out, ...
    var_vs_channel_x, var_vs_channel_y,sensor_force,sensor_force_arrival,sensor_force_peak,sensor_force_peak_time, ...
    sensor_vfrac,sensor_stressyy,wave_measurement,wave_time,wave_fit,wave_fit_time, ...
    init_overlap_mean_vs_depth, sensor_stress_peak_time, sensor_stress_peak_val,sensor_stress_yz, ...
    store_max_height_pdiff,height_out, shockspeedslope_exp, init_total_overlap_mean_vs_depth,store_heights_bin,radial_sensor_pos, ...
    sensor_overlap, sensor_overlap_peak, sensor_overlap_peak_time,Fmag_shocktrack_maxoverlap,Fmag_shocktrack_maxoverlap_time, ...
    Fmag_shocktrack_maxvelocity,Fmag_shocktrack_maxvelocity_time, F_overlaptrack_mean, F_overlaptrack_mean_time, F_velocitytrack_mean,vf_pdiff_voro, ...
    vf_rowavg_final, vf_rowstd_final, vf_colavg_final, vf_colstd_final,wspacing, F_overlaptrack_mean_vsdepth, F_velocitytrack_mean_vsdepth,F_forcetrack_mean_vsdepth, ...
    F_mag_avg_time,F_mag_avg,loft_depth_track] ...
    = HarvestAndPlotFun_production_loadonly(mat_file,plotpathdir,evalval,comptype,workingdir,vel,tot_deltas_onoff,SF,r_p)

% ************************************************************************
% This function performs the averaging procedure on the already harvested *.mat file
% Many of the outputs are vestigial and/or experimental so I will just describe the most
% relevant outputs here. This code would look way better (and probably run faster)
% if it were object oriented.
%
% % ----------- output ----------- %
% - channel - radial position of the sensors
% - depthplot - vertical position of the sensors
% - col/row avg/std - either the average or the standard deviation of sensor data,
% with respect to the columns (averaged across all radial positions at a single
% depth) or the rows (averaged across all the depths at a single radial position)
% - vx, vy, vz - particle velocities according to direction
% - Fx, Fy, Fz, Fmag - particle forces in various directions and the magnitude
% - timeplot - time array (s)
% - vf_rowavg_init - packing fraction at each depth (manually computed)
% - vf_pdiff_voro - packing fraction at eacth depth (via voronoi tesselation)
% - F_shocktrack_mean - average max particle force in wavefront vs channel position
% - init_overlap_mean_vs_depth - initial average particle overlap vs depth
% - F_overlaptrack_mean - average maximum particle overlap in the wavefront at each
% radial position
% - F_overlaptrack_mean_time - avg. max. particle overlap in the wavefront vs time
% - F_velocitytrack_mean - avg. max. particle velocity in the wavefront vs radial pos.
% - F_velocitytrack_mean_vsdepth - avg. max particle overlap in the wavefront at a
% single location
% - F_forcetrack_mean_vsdepth - avg. max particle force in the wavefront at a single
% location
%
% % ----------- input ----------- %
% - mat_file - harvested *.mat file containing virtual sensor averaged state info vs time
% - plotpathdir - directory to save plots in
% - evalval - label for this specific test within the cases considered
% - comptype - parameter considered for these cases
% - workingdir - directory with the code
% - vel - piston velocity (m/s)
% - tot_deltas_onoff - vestigial flag
% - SF - scaling flag
% - r_p - particle radius (m)
%
% ************************************************************************




% find the grain line
%[grainline, channellength] = findGrainLineEvery(filename,0);


% time and date information and strings
datetime.setDefaultFormats('default','yyyy-MM-dd')
currenttime = datetime;
formatOut = 'mmddyyyy';
timestring = datestr(currenttime,formatOut);
% grainstring = num2str(round(100*grainline));
% channelstring = num2str(channellength);
% typestring = num2str(evalval);


%% load existing data or harvest it

% see if any files exist that match the naming convention
%qfiles = dir(fullfile(harvestpath, storestringshort));

% if the file exists, just load that data and update the store string
% if ~isempty(qfiles)
%     % if there is more than 1 file, pause
%     checkstructsize = size(qfiles);
%     if checkstructsize(1) > 1
%         keyboard
%     else
%         storestring = qfiles.name;
load(mat_file);



plotpathdir = [workingdir 'ProductionPlots/' comptype '/'];


% redefine the sensor cell/depth here (it gets saved from the first run,
% and then uses that if not updated here
% sz = 3;%10;
% sr = 20;%82;%20;
% approximate sensor to be in the same location (1m length, 10 cm depth)

if SF
    sz = floor(SF*0.01/wspacing);%10cm;
    sr = floor(SF*0.01/wspacing);%82;%20;
else
    sz = floor(0.1/wspacing);%10;
    sr = floor(1/wspacing);%82;%20;
end
% set the radial sensor indicator
radial_sensor_pos = sr*wspacing;

% make sure plot path is current
plotpath = [plotpathdir num2str(evalval) '/'];
plotsonind = 1;

% define region bounds (done from plot comparison and analysis)
% wall region - portion down channel affected by wall
%[(end wall - wall reg), wall]
wallbuff = 1.50; % m
floorbuff = .18; % m
wallregstart = 100*(ystop - wallbuff); % cm
floorregstart = 100*(zstart + floorbuff); % cm

% define length down channel for plotting against
channel = 100 * (windowrange(1):wspacing:(windowrange(end)-wspacing));
depthplot = 100 * (zstart:wspacing:(zstop-wspacing));
% define depth and channel length bounds
depthstart = 1;
depthstop = znumgrids;


% define the regions in terms of grid numbers
regionstopradial = find(channel >= wallregstart,1,'first') - 1;
% regionstopdepth = find(depthplot >= floorregstart,1,'first') - 1;

range_prop = 55;
range_exp = 30;


grainstring = num2str(round(100*grainline));
channelstring = num2str(channellength);
porstring = num2str(round(100*porosity));

%% post process loop
% find time index corresponding to the last 0.5 seconds of sim time
last_pointfive_ind = find(timeplot > (timeplot(end)-0.5),1,'first');

%% initial conditions

% find the initial delta vs depth
% ignore the walls
init_overlaps_vs_depth = cell(znumgrids,1);
init_overlap_mean_vs_depth = zeros(znumgrids,1);
init_overlap_std_vs_depth= zeros(znumgrids,1);

init_total_overlaps_vs_depth = cell(znumgrids,1);
init_total_overlap_mean_vs_depth = zeros(znumgrids,1);
init_total_overlap_std_vs_depth= zeros(znumgrids,1);

init_sensor_overlaps_vs_depth = cell(znumgrids,1);
init_sensor_overlap_mean_vs_depth = zeros(znumgrids,1);
init_sensor_overlap_std_vs_depth= zeros(znumgrids,1);


for jj = 1:znumgrids
    for ii = 10:(ynumgrids-10) % wall exclusion
        init_overlaps_vs_depth{jj,1} = [init_overlaps_vs_depth{jj,1}, store_delta{jj,ii}(1)];
        if tot_deltas_onoff
            init_total_overlaps_vs_depth{jj,1} = [init_total_overlaps_vs_depth{jj,1}, store_delta_total{jj,ii}(1)];
        end
            init_sensor_overlaps_vs_depth{jj,1} = [init_sensor_overlaps_vs_depth{jj,1}, store_psensor_delta{jj,1}(1)];

    end


    % and then take the average
    % nan protection
    temp_overlap_nans = ~isnan(init_overlaps_vs_depth{jj,1});
    init_overlap_mean_vs_depth(jj,1) = mean(init_overlaps_vs_depth{jj,1}(temp_overlap_nans));
    init_overlap_std_vs_depth(jj,1) = std(init_overlaps_vs_depth{jj,1}(temp_overlap_nans));



    if tot_deltas_onoff
        temp_total_overlap_nans = ~isnan(init_total_overlaps_vs_depth{jj,1});
        init_total_overlap_mean_vs_depth(jj,1) = mean(init_total_overlaps_vs_depth{jj,1}(temp_total_overlap_nans));
        init_total_overlap_std_vs_depth(jj,1) = std(init_total_overlaps_vs_depth{jj,1}(temp_total_overlap_nans));
    end

    init_sensor_overlap_mean_vs_depth(jj,1) = mean(init_sensor_overlaps_vs_depth{jj,1});
    init_sensor_overlap_std_vs_depth(jj,1) = std(init_sensor_overlaps_vs_depth{jj,1});



end

% force
[Fmag_rowavg_init, Fmag_rowstd_init, Fmag_colavg_init, Fmag_colstd_init] = ...
    findInitialConditions(store_Fmag,ynumgrids,znumgrids);
[Fx_rowavg_init, Fx_rowstd_init, Fx_colavg_init, Fx_colstd_init] = ...
    findInitialConditions(store_Fx,ynumgrids,znumgrids);
[Fy_rowavg_init, Fy_rowstd_init, Fy_colavg_init, Fy_colstd_init] = ...
    findInitialConditions(store_Fy,ynumgrids,znumgrids);
[Fz_rowavg_init, Fz_rowstd_init, Fz_colavg_init, Fz_colstd_init] = ...
    findInitialConditions(store_Fz,ynumgrids,znumgrids);


% velocity
[vmag_rowavg_init, vmag_rowstd_init, vmag_colavg_init, vmag_colstd_init] = ...
    findInitialConditions(store_Vmag,ynumgrids,znumgrids);
[vy_rowavg_init, vy_rowstd_init, vy_colavg_init, vy_colstd_init] = ...
    findInitialConditions(store_VY,ynumgrids,znumgrids);
[vx_rowavg_init, vx_rowstd_init, vx_colavg_init, vx_colstd_init] = ...
    findInitialConditions(store_VX,ynumgrids,znumgrids);
[vz_rowavg_init, vz_rowstd_init, vz_colavg_init, vz_colstd_init] = ...
    findInitialConditions(store_VZ,ynumgrids,znumgrids);


% volume fraction
% [vf_rowavg_init, vf_rowstd_init, vf_colavg_init, vf_colstd_init] = ...
%     findInitialConditions_vf(store_vf_mean_post,ynumgrids,znumgrids);

[vf_rowavg_init, vf_rowstd_init, vf_colavg_init, vf_colstd_init] = ...
    findInitialConditions_vf(store_voro_bins,ynumgrids,znumgrids);


[z_end,y_end] = size(voro_bins_end);
[vf_rowavg_final, vf_rowstd_final, vf_colavg_final, vf_colstd_final] = ...
    findFinalConditions_vf(voro_bins_end,y_end,z_end);
% out put these and make a plot. also need depthplot/axis descriptor



% just take the output that is the same length as the initial depthplot so
% we can use the same array for plotting
% the above works well if the bed height stays the same or increases, but
% not if the bed is compacted. append NaN to end of array to make equal
% length as original depthplot if that happens
if z_end < znumgrids
    vf_rowavg_final = [vf_rowavg_final, nan(1,(znumgrids - z_end))];
    vf_rowstd_final = [vf_rowstd_final, nan(1,(znumgrids - z_end))];
    vf_colavg_final = [vf_colavg_final, nan(1,(znumgrids - z_end))];
    vf_colstd_final = [vf_colstd_final, nan(1,(znumgrids - z_end))];
else
    vf_rowavg_final = vf_rowavg_final(1:znumgrids);
    vf_rowstd_final = vf_rowstd_final(1:znumgrids);
    vf_colavg_final = vf_colavg_final(1:znumgrids);
    vf_colstd_final = vf_colstd_final(1:znumgrids);
end

vf_pdiff_voro = 100*(vf_rowavg_final - vf_rowavg_init) ./ vf_rowavg_init;

% zeros array for plotting depth error bars
zzz = zeros(1,znumgrids);


%% plot the initial conditions
channelstring = num2str(channellength);
if 0
% column averages
h_init_col = figure;
subplot(1,3,1)
errorbar(channel,Fmag_colavg_init,Fmag_colstd_init,'--o','LineWidth',2);
hold on
grid on
errorbar(channel,Fx_colavg_init,Fx_colstd_init,'--o','LineWidth',2);
errorbar(channel,Fy_colavg_init,Fy_colstd_init,'--o','LineWidth',2);
errorbar(channel,Fz_colavg_init,Fz_colstd_init,'--o','LineWidth',2);
xlabel('Radial Position [cm]')
%set(gca, 'YScale', 'log')
ylabel('Force [N] (log)')
legend('Magnitude','X','Y','Z')
title('Initial force vs radial position')
set(gca,'fontsize',14)

%figure
subplot(1,3,2)
errorbar(channel,100*vmag_colavg_init,100*vmag_colstd_init,'--o','LineWidth',2);
hold on
grid on
errorbar(channel,100*vx_colavg_init,100*vx_colstd_init,'--o','LineWidth',2);
errorbar(channel,100*vy_colavg_init,100*vy_colstd_init,'--o','LineWidth',2);
errorbar(channel,100*vz_colavg_init,100*vz_colstd_init,'--o','LineWidth',2);
xlabel('Radial Position [cm]')
ylabel('Velocity [cm/s]')
legend('Magnitude','X','Y','Z')
title('Initial velocity vs radial position')
set(gca,'fontsize',14)

%figure
subplot(1,3,3)
errorbar(channel,100*vf_colavg_init,100*vf_colstd_init,'--o','LineWidth',2);
hold on
grid on
xlabel('Radial Position [cm]')
ylabel('Volume Fraction [%]')
title('Initial volume fraction vs radial position')
set(gca,'fontsize',14)

sgtitle(['Depth averaged initial condition plots - ' channelstring 'm channel'])

set(h_init_col, 'Units', 'pixels', 'Position', [100, 300, 1400, 400]);


% row averages
h_init_row = figure;
subplot(1,3,1)
errorbar(Fmag_rowavg_init,depthplot,zzz,zzz,Fmag_rowstd_init,Fmag_rowstd_init,'--o','LineWidth',2);
hold on
grid on
errorbar(Fx_rowavg_init,depthplot,zzz,zzz,Fx_rowstd_init,Fx_rowstd_init,'--o','LineWidth',2);
errorbar(Fy_rowavg_init,depthplot,zzz,zzz,Fy_rowstd_init,Fy_rowstd_init,'--o','LineWidth',2);
errorbar(Fz_rowstd_init,depthplot,zzz,zzz,Fz_rowavg_init,Fz_rowavg_init,'--o','LineWidth',2);
ylabel('Depth [cm]')
set(gca, 'XScale', 'log')
set(gca,'YDir','reverse')
xlabel('Force [N] (log)')
legend('Magnitude','X','Y','Z')
title('Depth vs Initial force')
set(gca,'fontsize',14)

%figure
subplot(1,3,2)
errorbar(100*vmag_rowavg_init,depthplot,100*vmag_rowstd_init,'horizontal','--o','LineWidth',2);
hold on
grid on
errorbar(100*vx_rowavg_init,depthplot,zzz,zzz,100*vx_rowstd_init,100*vx_rowstd_init,'--o','LineWidth',2);
errorbar(100*vy_rowavg_init,depthplot,zzz,zzz,100*vy_rowstd_init,100*vy_rowstd_init,'--o','LineWidth',2);
errorbar(100*vz_rowavg_init,depthplot,zzz,zzz,100*vz_rowstd_init,100*vz_rowstd_init,'--o','LineWidth',2);
ylabel('Depth [cm]')
xlabel('Velocity [cm/s]')
set(gca,'YDir','reverse')
legend('Magnitude','X','Y','Z')
title('Depth vs Initial velocity')
set(gca,'fontsize',14)

%figure
subplot(1,3,3)
errorbar(100*vf_rowavg_init,depthplot,zzz,zzz,100*vf_rowstd_init,100*vf_rowstd_init,'--o','LineWidth',2);
hold on
grid on
ylabel('Depth [cm]')
xlabel('Volume Fraction [%]')
set(gca,'YDir','reverse')
title('Depth vs Initial volume fraction')
set(gca,'fontsize',14)

sgtitle(['Radially averaged initial condition plots - '  channelstring 'm channel - ' grainstring ' cm fill' ])

set(h_init_row, 'Units', 'pixels', 'Position', [100, 300, 1400, 400]);

end

%% define analysis region based off of initial conditions

%use the average initial force to find the depth cutoff

%regionstopdepth = find(depthplot >= floorregstart,1,'first') - 1;
regionstartdepth = 1; % exculde top layer

% for the pressurized case, need to remove the upper layers that get
% corrupted
if evalval == "20 m regolith" || evalval == "1 m regolith" || evalval ==  "5 m regolith" || ...
        evalval == "10 m regolith"
    regionstartdepth = 9;

end

% stop 5 particle diameters before the all (1.25 cm)
regionstopdepth = znumgrids - 3;

if comptype == "decay_boulder_layers"
    regionstopdepth = znumgrids;
end

%% first row dilation inspection

% subplot with
% plot of dilation vs time for each channel position
% plot of total dilated vloume vs time

% loop for creating dilated volume
% dilatedbinsvstime = zeros(1,Ltp)
% for jj = 1:ynumgrids
% end

% just evaluate dilation vs time for each channel position in the top row
% first
% %keyboard
% colcolorfull = jet(ynumgrids);
% % rowcolorfull = jet(length(1:13));
% figure
% axes('NextPlot','add')
%
% for jj = 1:ynumgrids
%
%     channelstr = ['Y = ' num2str(channel(jj)) ' cm'];
%     yyaxis left
%     hold on
%     if ssdilationind_full(1,jj) == 50
%         subplot(1,2,1)
%         plot(timeplot,store_vf_pdiff{1,jj},'--o','DisplayName',[channelstr],'Color',colcolorfull(jj,:))
%         grid on
%         ylabel('Percent Dilation')
%         xlabel('Time [s]')
%         title('Dilation over time for each radial position (full)')
%         legend('show','Location','Best')
%     elseif ssdilationind_full(1,jj) == -50
%         subplot(1,2,2)
%         plot(timeplot,store_vf_pdiff{1,jj},'--o','DisplayName',[channelstr],'Color',colcolorfull(jj,:))
%         grid on
%         ylabel('Percent Dilation')
%         xlabel('Time [s]')
%         title('Dilation over time for each radial position (Partial)')
%         legend('show','Location','Best')
%     end
%     grid on
%     ylabel('Percent Dilation')
%     xlabel('Time [s]')
%     %set(gca, 'YScale', 'log')
%     %title('Dilation over time for each radial position')
%     %set(gca,'fontsize',14)
%
% end
%
% sgtitle('Dilation over time for the top row')


%% steady state averaging by row and column, for export in combined plots
% row and column averaging for the different groupings

% % cut off into regions - a for analysis region
% ssdilationind_full_a = ssdilationind_full(2:regionstopdepth,1:regionstopradial);
% dilationpdiff_SS_a = dilationpdiff_SS(2:regionstopdepth,1:regionstopradial);
% depthplot_a = depthplot(regionstartdepth:regionstopdepth);
% channel_a = channel(1:regionstopradial);
%
%
% % --- dilation --- %
% % find the position of all cells which have dilated
% [dilrowind,dilcolind] = find(ssdilationind_full_a == 50);
% dilrowunq = unique(dilrowind); % get the rows that have any dilation
% dilcolunq = unique(dilcolind); % get the columns that have any dilation
% dilrows = length(dilrowunq);
% dilcols = length(dilcolunq);
% % storage
% ss_dilation_col_avg = zeros(1, dilcols);
% ss_dilation_col_std = zeros(1, dilcols);
% ss_dilation_col_avg_yval = zeros(1, dilcols);
% ss_dilation_row_avg = zeros(1, dilrows);
% ss_dilation_row_std = zeros(1, dilrows);
% ss_dilation_row_avg_yval = zeros(1, dilrows);
% for j = 1:dilrows
%
%     tempdil = zeros(1,(length(dilrowind(dilrowind == dilrowunq(j)))));
%     tempctr = 1;
%
%     for k = 1:length(dilcolind)
%
%
%
%         if dilrowind(k) == dilrowunq(j)
%             tempdil(1,tempctr) =  dilationpdiff_SS_a(dilrowind(k),dilcolind(k));
%             tempctr = tempctr + 1;
%
%         end
%     end
%
%     % if there are at least 2 cells that dilated, include them for
%     % averaging
%     if tempctr > 2
%         ss_dilation_row_avg(1,j) = nanmean(tempdil);
%         ss_dilation_row_std(1,j) = nanstd(tempdil);
%         ss_dilation_row_avg_yval(1,j) = depthplot_a(dilrowunq(j));
%     else
%         ss_dilation_row_avg(1,j) = NaN;
%         ss_dilation_row_std(1,j) = NaN;
%         ss_dilation_row_avg_yval(1,j) = NaN;
%     end
%
% end
%
% for j = 1:dilcols
%
%     tempdil = zeros(1,(length(dilcolind(dilcolind == dilcolunq(j)))));
%     tempctr = 1;
%
%     for k = 1:length(dilcolind)
%
%
%
%         if dilcolind(k) == dilcolunq(j)
%             tempdil(1,tempctr) =  dilationpdiff_SS_a(dilrowind(k),dilcolind(k));
%             tempctr = tempctr + 1;
%
%         end
%     end
%
%
%     if tempctr > 2
%         ss_dilation_col_avg(1,j) = nanmean(tempdil);
%         ss_dilation_col_std(1,j) = nanstd(tempdil);
%         ss_dilation_col_avg_yval(1,j) = channel_a(dilcolunq(j));
%     else
%         ss_dilation_col_avg(1,j) = NaN;
%         ss_dilation_col_std(1,j) = NaN;
%         ss_dilation_col_avg_yval(1,j) = NaN;
%     end
%
%
% end
%
%
% % --- partial --- %
% [partrowind,partcolind] = find(ssdilationind_full_a == -50);
% partrowunq = unique(partrowind);
% partcolunq = unique(partcolind);
% partrows = length(partrowunq);
% partcols = length(partcolunq);
% ss_part_col_avg = zeros(1, partcols);
% ss_part_col_std = zeros(1, partcols);
% ss_part_col_avg_yval = zeros(1, partcols);
% ss_part_row_avg = zeros(1, partrows);
% ss_part_row_std = zeros(1, partrows);
% ss_part_row_avg_yval = zeros(1, partrows);
% for j = 1:partrows
%
%     temppart = zeros(1,(length(partrowind(partrowind == partrowunq(j)))));
%     tempctr = 1;
%
%     for k = 1:length(partcolind)
%
%
%
%         if partrowind(k) == partrowunq(j)
%             temppart(1,tempctr) =  dilationpdiff_SS_a(partrowind(k),partcolind(k));
%             tempctr = tempctr + 1;
%
%         end
%     end
%
%     ss_part_row_avg(1,j) = nanmean(temppart);
%     ss_part_row_std(1,j) = nanstd(temppart);
%     ss_part_row_avg_yval(1,j) = depthplot_a(partrowunq(j));
%
% end
%
% for j = 1:partcols
%
%     temppart = zeros(1,(length(partcolind(partcolind == partcolunq(j)))));
%     tempctr = 1;
%
%     for k = 1:length(partcolind)
%
%
%
%         if partcolind(k) == partcolunq(j)
%             temppart(1,tempctr) =  dilationpdiff_SS_a(partrowind(k),partcolind(k));
%             tempctr = tempctr + 1;
%
%         end
%     end
%
%     ss_part_col_avg(1,j) = nanmean(temppart);
%     ss_part_col_std(1,j) = nanstd(temppart);
%     ss_part_col_avg_yval(1,j) = channel_a(partcolunq(j));
%
%
% end
%
%
% % --- compaction --- %
% [compactrowind,compactcolind] = find(ssdilationind_full == -100);
% compactrowunq = unique(compactrowind);
% compactcolunq = unique(compactcolind);
% compactrows = length(compactrowunq);
% compactcols = length(compactcolunq);
% ss_compact_col_avg = zeros(1, compactcols);
% ss_compact_col_std = zeros(1, compactcols);
% ss_compact_col_avg_yval = zeros(1, compactcols);
% ss_compact_row_avg = zeros(1, compactrows);
% ss_compact_row_std = zeros(1, compactrows);
% ss_compact_row_avg_yval = zeros(1, compactrows);

% for j = 1:dilrows
%
%     tempcompact = zeros(1,(length(compactrowind(compactrowind == compactrowunq(j)))));
%     tempctr = 1;
%
%
%
%     for k = 1:length(unique(dilcolind))
%
%
%
%         if compactrowind(k) == compactrowunq(j)
%             tempcompact(1,tempctr) =  dilationpdiff_SS(compactrowind(k),compactcolind(k));
%             tempctr = tempctr + 1;
%
%         end
%     end
%
%     ss_compact_row_avg(1,j) = nanmean(tempcompact);
%     ss_compact_row_std(1,j) = nanstd(tempcompact);
%     ss_compact_row_avg_yval(1,j) = depthplot(compactrowunq(j));
%
% end
%
% for j = 1:dilcols
%
%     tempcompact = zeros(1,(length(compactcolind(dilcolind == compactcolunq(j)))));
%     tempctr = 1;
%
%     for k = 1:length(unique(compactcolind))
%
%
%
%         if compactcolind(k) == compactcolunq(j)
%             tempcompact(1,tempctr) =  dilationpdiff_SS(compactrowind(k),compactcolind(k));
%             tempctr = tempctr + 1;
%
%         end
%     end
%
%     ss_compact_col_avg(1,j) = nanmean(tempcompact);
%     ss_compact_col_std(1,j) = nanstd(tempcompact);
%     ss_compact_col_avg_yval(1,j) = channel(compactcolunq(j));
%
%
% end
%
% zdil = zeros(1,dilrows);

% if plotsonind
% figure
% subplot(1,2,1)
% errorbar(ss_dilation_col_avg_yval,ss_dilation_col_avg,ss_dilation_col_std,'--o','LineWidth',2)
% hold on
% grid on
% errorbar(ss_part_col_avg_yval,ss_part_col_avg,ss_part_col_std,'--o','LineWidth',2)
% %errorbar(ss_compact_col_avg_yval,ss_compact_col_avg,ss_compact_col_std,'--o','LineWidth',2)
% xlabel('Channel Position [cm]')
% ylabel('Percent Dilation')
% legend('Dilation','Partial Dilation','Location','Best')
% title('Depth averaged steady state dilation')
% set(gca,'fontsize',14)
%
% subplot(1,2,2)
% errorbar(ss_dilation_row_avg,ss_dilation_row_avg_yval,zdil,zdil,ss_dilation_row_std,ss_dilation_row_std,'--o','LineWidth',2)
% hold on
% grid on
% %errorbar(ss_part_row_avg,ss_part_row_avg_yval,zzz,zzz,ss_part_row_std,ss_part_row_std,'--o','LineWidth',2)
% %errorbar(ss_compact_row_avg_yval,ss_compact_row_avg,ss_compact_row_std,'--o','LineWidth',2)
% ylabel('Depth [cm]')
% xlabel('Percent Dilation')
% set(gca,'YDir','reverse')
% legend('Dilation','Partial Dilation','Location','Best')
% title('Depth vs steady state dilation')
% set(gca,'fontsize',14)
%
%
% sgtitle(['Steady State Dilations - ' channelstring 'm channel - ' grainstring ' cm fill' ])
%
% end
%% old
% stuff i haven't deleted yet but is currently in output so keeping for now
% max average particle velocity encountered at each depth
maxofmaxvel = (max(vy_max.')).';
percentofmaxvel = (vy_max - maxofmaxvel)./vy_max;


maxofmaxF = (max(Fmag_max.')).';
percentofmaxF = (Fmag_max - maxofmaxF)./Fmag_max;


%% find the dilation bins over time
% analysis region
% [dilationpdiff_SS_tp_f_a,dilationpdiff_SS_tp_p_a] = ...
%     findBinPlots(store_vf_pdiff,1,regionstopradial,1,regionstopdepth, ...
%     depthplot,'Analysis',Ltp,timeplot,channelstring,grainstring);
% % floor region
% [dilationpdiff_SS_tp_f_f,dilationpdiff_SS_tp_p_f] = ...
%     findBinPlots(store_vf_pdiff,1,regionstopradial,regionstopdepth,znumgrids, ...
%     depthplot,'Floor',Ltp,timeplot,channelstring,grainstring);
% % wall region
% [dilationpdiff_SS_tp_f_w,dilationpdiff_SS_tp_p_w] = ...
%     findBinPlots(store_vf_pdiff,regionstopradial,ynumgrids,1,regionstopdepth, ...
%     depthplot,'Wall',Ltp,timeplot,channelstring,grainstring);
% %keyboard

%% dilation over time - depth
%
% % uses store_vf_pdiff
% % make a nested loop that goes through each depth
%     % and make an eveerage dilation vs time at this depth across all the
%     % radial positions
%
% vfrac_vs_time_mean_depth = zeros(znumgrids,Ltp);
% vfrac_vs_time_std_depth = zeros(znumgrids,Ltp);
%
%
% for kk = 1:Ltp
%
%     temp_vf_vs_time_depth = zeros(1,ynumgrids);
%
%     for jj = 1:znumgrids
%
%         for ii = 1:ynumgrids
%
%         temp_vf_vs_time_depth(1,ii) = store_vf_pdiff{jj,ii}(kk);
%
%         end % end y grids
%
%         vfrac_vs_time_mean_depth(jj,kk) = nanmean(temp_vf_vs_time_depth);
%         vfrac_vs_time_std_depth(jj,kk) = nanstd(temp_vf_vs_time_depth);
%
%     end % end z grids
%     %keyboard
% end % end times
%
% h_vf_vs_time_depth = figure;
% axes('NextPlot','add')
% plotcolor = turbo(znumgrids);
%
% for jj = 2:znumgrids
%
%     subplot(1,1,1)
%     hold on
%     plot(timeplot,(vfrac_vs_time_mean_depth(jj,:)),'.','DisplayName',['Depth ' num2str(depthplot(jj)) ' cm'],'Color',plotcolor(jj,:))
%     grid on
%     xlabel('time [s]')
%     ylabel('Magnditude of average percent change in volume fraction [%]')
%
% end
% legend('show','location','best')
% title(['V_{im} = ' num2str(evalval), ' [m/s]'])
%
% %keyboard

%% dilation over time - radial

% uses store_vf_pdiff
% make a nested loop that goes through each depth
    % and make an eveerage dilation vs time at this depth across all the
    % radial positions

vfrac_vs_time_mean_rad = zeros(ynumgrids,Ltp);
vfrac_vs_time_std_rad = zeros(ynumgrids,Ltp);


% for kk = 1:Ltp
%
%     temp_vf_vs_time_rad = zeros(1,znumgrids);
%
%     % remove the top and bottom row from consideration
%     for jj = 1:ynumgrids
%
%         for ii = 1:znumgrids
%
%         temp_vf_vs_time_rad(1,ii) = store_vf_pdiff{ii,jj}(kk);
%
%         end % end z grids
%
%         vfrac_vs_time_mean_rad(jj,kk) = nanmean(temp_vf_vs_time_rad);%(2:(znumgrids-1)));
%         vfrac_vs_time_std_rad(jj,kk) = nanstd(temp_vf_vs_time_rad);%(2:(znumgrids-1)));
%
%     end % end y grids
%     %keyboard
% end % end times


plotcolor = turbo(ynumgrids);


% h_end_variance_vs_channel = figure;
% axes('NextPlot','add')


var_vs_channel_x = zeros(1,length(2:ynumgrids-5));
var_vs_channel_y = zeros(1,length(2:ynumgrids-5));
% initialize for storage
sensor_vfrac = 0;

%neglect the first radial index since the wave is initiated here
% for jj = 2:ynumgrids-5
%
%     % get derivative of the high res stuff, want to know chang at end time
%     vf_lowres = vfrac_vs_time_mean_rad(jj,S_lowres:Ltp);
%     %time_lowres = timeplot(S_lowres:Ltp);
%
%     %dtdvf = finite_diff_d1_fun(vf_lowres,timestep2*simtimestep);
%
%     if  vfrac_vs_time_mean_rad(jj,Ltp) >= 0    %jj < 30
%
%
%
%     else
%
%
%
%
%     end
%
%         last_volfracs = zeros(1,length(4:znumgrids-1));
%         deltaphis = zeros(1,length(4:znumgrids-1));
%
%         % neglect the floor since frozen particles are there
%     for kk = 4:(znumgrids-1)
%         temp_last_vfs = store_vf_mean_post{kk,jj}((end-20):end);
%         last_volfracs(1,kk-3) = std(temp_last_vfs);
%
%         % find the delta phi as defined in mathias for the last time steps,
%         % then put their average here
%         deltaphis(1,kk-3) = mean(abs(temp_last_vfs - mean(temp_last_vfs)));
%
%         % log for sensor output
%         if kk == sz &&  jj == sr
%             sensor_vfrac = store_vf_mean_post{kk,jj};
%         end
%
%     end
%
%
%         set(0, 'CurrentFigure', h_end_variance_vs_channel)
%         subplot(1,1,1)
%         hold on
%         plot(channel(jj)/100,10^3*mean(last_volfracs),'ko')
%         %semilogy(channel(jj)/100,abs(mean(vfrac_vs_time_mean_rad(jj,(end-10):end))),'ko')
%         grid on
%         xlabel('R [m]')
%         ylabel('10^3\sigma_{\phi}')
%         %title(['R < 30 cm, V_{im} = ' num2str(evalval), ' [m/s]'])
%         title(['Average \sigma\phi_{ss}, V_{im} = ' num2str(evalval), ' [m/s]'])
%         set(gca,'fontsize',14)
%         set(gca, 'YScale', 'log')
%
%
%         var_vs_channel_x(1,jj-1) =channel(jj)/100;
%         var_vs_channel_y(1,jj-1) = 10^3*mean(last_volfracs);
%
%
%
% %     keyboard
% end


%keyboard
%% row and column averages for varying types (full dilation, partial, compaction)

% % magnitude
% [Fmag_rowavg,Fmag_colavg,vmag_rowavg,vmag_colavg,vfpdiff_rowavg, vfpdiff_colavg, ...
% Fmag_rowstd,Fmag_colstd,vmag_rowstd,vmag_colstd,vfpdiff_rowstd,vfpdiff_colstd, ...
% Fmag_rowavg_partial,Fmag_colavg_partial,vmag_rowavg_partial,vmag_colavg_partial,vfpdiff_rowavg_partial,vfpdiff_colavg_partial, ...
% Fmag_rowstd_partial,Fmag_colstd_partial,vmag_rowstd_partial,vmag_colstd_partial,vfpdiff_rowstd_partial,vfpdiff_colstd_partial, ...
% Fmag_rowavg_compact,Fmag_colavg_compact,vmag_rowavg_compact,vmag_colavg_compact,vfpdiff_rowavg_compact,vfpdiff_colavg_compact, ...
% Fmag_rowstd_compact,Fmag_colstd_compact,vmag_rowstd_compact,vmag_colstd_compact,vfpdiff_rowstd_compact,vfpdiff_colstd_compact] = ...
%     findWindowAverages(ynumgrids,znumgrids,Ltp,store_vf_pdiff,store_Vmag,store_Fmag,ssdilationind_full);
% 
% % y
% [Fy_rowavg,Fy_colavg,vy_rowavg,vy_colavg,vfpdiff_rowavg, vfpdiff_colavg, ...
% Fy_rowstd,Fy_colstd,vy_rowstd,vy_colstd,vfpdiff_rowstd,vfpdiff_colstd, ...
% Fy_rowavg_partial,Fy_colavg_partial,vy_rowavg_partial,vy_colavg_partial,vfpdiff_rowavg_partial,vfpdiff_colavg_partial, ...
% Fy_rowstd_partial,Fy_colstd_partial,vy_rowstd_partial,vy_colstd_partial,vfpdiff_rowstd_partial,vfpdiff_colstd_partial, ...
% Fy_rowavg_compact,Fy_colavg_compact,vy_rowavg_compact,vy_colavg_compact,vfpdiff_rowavg_compact,vfpdiff_colavg_compact, ...
% Fy_rowstd_compact,Fy_colstd_compact,vy_rowstd_compact,vy_colstd_compact,vfpdiff_rowstd_compact,vfpdiff_colstd_compact] = ...
%     findWindowAverages(ynumgrids,znumgrids,Ltp,store_vf_pdiff,store_VY,store_Fy,ssdilationind_full);
% 
% % x
% [Fx_rowavg,Fx_colavg,vx_rowavg,vx_colavg,vfpdiff_rowavg, vfpdiff_colavg, ...
% Fx_rowstd,Fx_colstd,vx_rowstd,vx_colstd,vfpdiff_rowstd,vfpdiff_colstd, ...
% Fx_rowavg_partial,Fx_colavg_partial,vx_rowavg_partial,vx_colavg_partial,vfpdiff_rowavg_partial,vfpdiff_colavg_partial, ...
% Fx_rowstd_partial,Fx_colstd_partial,vx_rowstd_partial,vx_colstd_partial,vfpdiff_rowstd_partial,vfpdiff_colstd_partial, ...
% Fx_rowavg_compact,Fx_colavg_compact,vx_rowavg_compact,vx_colavg_compact,vfpdiff_rowavg_compact,vfpdiff_colavg_compact, ...
% Fx_rowstd_compact,Fx_colstd_compact,vx_rowstd_compact,vx_colstd_compact,vfpdiff_rowstd_compact,vfpdiff_colstd_compact] = ...
%     findWindowAverages(ynumgrids,znumgrids,Ltp,store_vf_pdiff,store_VX,store_Fx,ssdilationind_full);
% 
% % z
% [Fz_rowavg,Fz_colavg,vz_rowavg,vz_colavg,vfpdiff_rowavg, vfpdiff_colavg, ...
% Fz_rowstd,Fz_colstd,vz_rowstd,vz_colstd,vfpdiff_rowstd,vfpdiff_colstd, ...
% Fz_rowavg_partial,Fz_colavg_partial,vz_rowavg_partial,vz_colavg_partial,vfpdiff_rowavg_partial,vfpdiff_colavg_partial, ...
% Fz_rowstd_partial,Fz_colstd_partial,vz_rowstd_partial,vz_colstd_partial,vfpdiff_rowstd_partial,vfpdiff_colstd_partial, ...
% Fz_rowavg_compact,Fz_colavg_compact,vz_rowavg_compact,vz_colavg_compact,vfpdiff_rowavg_compact,vfpdiff_colavg_compact, ...
% Fz_rowstd_compact,Fz_colstd_compact,vz_rowstd_compact,vz_colstd_compact,vfpdiff_rowstd_compact,vfpdiff_colstd_compact] = ...
%     findWindowAverages(ynumgrids,znumgrids,Ltp,store_vf_pdiff,store_VZ,store_Fz,ssdilationind_full);

%% re do the averaging for velocity and force, based on all windows

% -- force -- %
% magnitude
[Fmag_rowavg_full,Fmag_colavg_full,Fmag_rowstd_full,Fmag_colstd_full] = ...
    findFullAverages(ynumgrids,znumgrids,Ltp,store_Fmag);

% xy magnitude
[Fmag_xy_rowavg_full,Fmag_xy_colavg_full,Fmag_xy_rowstd_full,Fmag_xy_colstd_full] = ...
    findFullAverages_xy(ynumgrids,znumgrids,Ltp,store_Fy,store_Fx);
% y
[Fy_rowavg_full,Fy_colavg_full,Fy_rowstd_full,Fy_colstd_full] = ...
    findFullAverages(ynumgrids,znumgrids,Ltp,store_Fy);
% x
[Fx_rowavg_full,Fx_colavg_full,Fx_rowstd_full,Fx_colstd_full] = ...
    findFullAverages(ynumgrids,znumgrids,Ltp,store_Fx);
% z
[Fz_rowavg_full,Fz_colavg_full,Fz_rowstd_full,Fz_colstd_full] = ...
    findFullAverages(ynumgrids,znumgrids,Ltp,store_Fz);

% -- velocity -- %
% magnitude
[vmag_rowavg_full,vmag_colavg_full,vmag_rowstd_full,vmag_colstd_full] = ...
    findFullAverages(ynumgrids,znumgrids,Ltp,store_Vmag);
% y
[vy_rowavg_full,vy_colavg_full,vy_rowstd_full,vy_colstd_full] = ...
    findFullAverages(ynumgrids,znumgrids,Ltp,store_VY);
% x
[vx_rowavg_full,vx_colavg_full,vx_rowstd_full,vx_colstd_full] = ...
    findFullAverages(ynumgrids,znumgrids,Ltp,store_VX);
% z
[vz_rowavg_full,vz_colavg_full,vz_rowstd_full,vz_colstd_full] = ...
    findFullAverages(ynumgrids,znumgrids,Ltp,store_VZ);

% -- volume fraction -- %
% [vf_rowavg_full,vf_colavg_full,vf_rowstd_full,vf_colstd_full] = ...
%     findFullAverages(ynumgrids,znumgrids,Ltp,store_vf_pdiff);




%% compute wave tracking - all windows - using averages across all depths at a given radial position

% find the first minimums from force over time at each depth (duration of
% shock period) - use shock max as the second time setp (t ~= 0)

% max is the maximum of the
shockmax = max(Fmag_rowavg_full(:,1:3).').';
%shockmax = Fmag_rowavg_full(:,2);
shockmaxstd = Fmag_rowstd_full(:,2);
shockmaxtime = zeros(znumgrids,1);
shockmaxtime(:,1) = timeplot(2);
shockmin = zeros(znumgrids,1);
shockminstd = zeros(znumgrids,1);
shockmintime = zeros(znumgrids,1);


for jj = 1:znumgrids

    shock_depthmax = shockmax(jj);
    shock_depthprior = shock_depthmax;
    minflag = 0;
    shockcounter = 2;

    while minflag < 1

        tempshockmin = Fmag_rowavg_full(jj,shockcounter);

        if tempshockmin > shock_depthprior
            shockmin(jj,1) = tempshockmin;
            shockminstd(jj,1) = Fmag_rowstd_full(jj,shockcounter);
            shockmintime(jj,1) = timeplot(shockcounter);
            minflag = minflag + 1;
        end

        shock_depthprior = tempshockmin;
        shockcounter = shockcounter + 1;

        if shockcounter > Ltp
            minflag = minflag + 1;
        end

    end

end

% shockperiodmaxtime = max(shockmintime);
% shockperiodmaxtimeind = floor(shockperiodmaxtime/(simtimestep*timestep1))+1;
 %spmtind = shockperiodmaxtimeind; %  shorten name

spmtind = shockcounter;


%% find upper and lower grouping from velocity maximums over shock period
% [shockmaxvels, shockmaxvelind] = max(vy_rowavg_full(:,1:shockperiodmaxtimeind).');
%
% aveind = mean(shockmaxvelind);
% grouper = zeros(1,znumgrids);
% for jj = 1:znumgrids
%     if shockmaxvelind(jj) < aveind
%         grouper(1,jj) = 1;
%     end
% end
% Llower = nnz(grouper);


%% find maximums

% -- force -- %
% magnitude
[max_Fmag_radial_fullw, max_Fmagtime_radial_fullw, max_Fmagdist_radial_fullw, max_Fmagstd_radial_fullw , ...
    max_Fmag_depth_fullw, max_Fmagtime_depth_fullw, max_Fmagdepth_fullw] ...
     = findFullMax(ynumgrids,znumgrids,spmtind, wspacing, windowrange,timeplot, ...
     Fmag_colavg_full,Fmag_colstd_full,Fmag_rowavg_full);
% x
[max_Fxradial_fullw, max_Fxtime_radial_fullw, max_Fxdist_radial_fullw, max_Fxstd_radial_fullw , ...
    max_Fx_depth_fullw, max_Fxtime_depth_fullw, max_Fxdepth_fullw] ...
     = findFullMax(ynumgrids,znumgrids,spmtind, wspacing, windowrange,timeplot, ...
     Fx_colavg_full,Fx_colstd_full,Fx_rowavg_full);
% y
[max_Fyradial_fullw, max_Fytime_radial_fullw, max_Fydist_radial_fullw, max_Fystd_radial_fullw , ...
    max_Fy_depth_fullw, max_Fytime_depth_fullw, max_Fydepth_fullw] ...
     = findFullMax(ynumgrids,znumgrids,spmtind, wspacing, windowrange,timeplot, ...
     Fy_colavg_full,Fy_colstd_full,Fy_rowavg_full);
% z
[max_Fzradial_fullw, max_Fztime_radial_fullw, max_Fzdist_radial_fullw, max_Fzstd_radial_fullw , ...
    max_Fz_depth_fullw, max_Fztime_depth_fullw, max_Fzdepth_fullw] ...
     = findFullMax(ynumgrids,znumgrids,spmtind, wspacing, windowrange,timeplot, ...
     Fz_colavg_full,Fz_colstd_full,Fz_rowavg_full);

% -- velocity -- %
% magnitude
[max_vmag_radial_fullw, max_vmagtime_radial_fullw, max_vmagdist_radial_fullw, max_vmagstd_radial_fullw , ...
    max_vmag_depth_fullw, max_vmagtime_depth_fullw, max_vmagdepth_fullw] ...
     = findFullMax(ynumgrids,znumgrids,spmtind, wspacing, windowrange,timeplot, ...
     vmag_colavg_full,vmag_colstd_full,vmag_rowavg_full);

% x
[max_vx_radial_fullw, max_vxtime_radial_fullw, max_vxdist_radial_fullw, max_vxstd_radial_fullw , ...
    max_vx_depth_fullw, max_vxtime_depth_fullw, max_vxdepth_fullw] ...
     = findFullMax(ynumgrids,znumgrids,spmtind, wspacing, windowrange,timeplot, ...
     vx_colavg_full,vx_colstd_full,vx_rowavg_full);
% y
[max_vy_radial_fullw, max_vytime_radial_fullw, max_vydist_radial_fullw, max_vystd_radial_fullw , ...
    max_vy_depth_fullw, max_vytime_depth_fullw, max_vydepth_fullw] ...
     = findFullMax(ynumgrids,znumgrids,spmtind, wspacing, windowrange,timeplot, ...
     vy_colavg_full,vy_colstd_full,vy_rowavg_full);
% z
[max_vz_radial_fullw, max_vztime_radial_fullw, max_vzdist_radial_fullw, max_vzstd_radial_fullw , ...
    max_vz_depth_fullw, max_vztime_depth_fullw, max_vzdepth_fullw] ...
     = findFullMax(ynumgrids,znumgrids,spmtind, wspacing, windowrange,timeplot, ...
     vz_colavg_full,vz_colstd_full,vz_rowavg_full);
% -- volume fracdtion -- %
% [max_vf_radial_fullw, max_vftime_radial_fullw, max_vfdist_radial_fullw, max_vfstd_radial_fullw , ...
%     max_vf_depth_fullw, max_vftime_depth_fullw, max_vfdepth_fullw] ...
%      = findFullMax(ynumgrids,znumgrids,spmtind, wspacing, windowrange,timeplot, ...
%      vf_colavg_full,vf_colstd_full,vf_rowavg_full);
%


%% plot over time
lengthstart = 1;
lengthstop = ynumgrids;

% llength_dil = nnz(sum(ssdilationind_binary(:,lengthstart:lengthstop)));
% dlength_dil = nnz(sum(ssdilationind_binary(depthstart:depthstop,:).'));
% 
% llength_partial = nnz(sum(vfpdiff_colavg_partial.'));
% dlength_partial = nnz(sum(vfpdiff_rowavg_partial.'));
% 
% llength_compact = nnz(sum(vfpdiff_colavg_compact.'));
% dlength_compact = nnz(sum(vfpdiff_rowavg_compact.'));

% llength_dil = nnz(sum(ssdilationind_binary(:,lengthstart:lengthstop)));
% dlength_dil = nnz(sum(ssdilationind_binary(regionstartdepth:regionstopdepth,:).'));
%
% llength_partial = nnz(sum(vfpdiff_colavg_partial.'));
% dlength_partial = nnz(sum(vfpdiff_rowavg_partial.'));
%
% llength_compact = nnz(sum(vfpdiff_colavg_compact.'));
% dlength_compact = nnz(sum(vfpdiff_rowavg_compact.'));


% generate plots
% --- magnitude --- %
%dilation

% [h_vsdepth_forceanddil] = findWindowPlots(dlength_dil,llength_dil,'dilation', 'magnitude', windowrange, wspacing, timeplot, ...
%     vfpdiff_colavg, vfpdiff_colstd, vfpdiff_rowavg, vfpdiff_rowstd, ...
%     vmag_colavg,vmag_colstd, vmag_rowavg, vmag_rowstd, ...
%     Fmag_colavg, Fmag_colstd, Fmag_rowavg, Fmag_rowstd, channelstring, grainstring,ynumgrids,znumgrids,Ltp);

% save figure
% saveas(h_dilationovertime, [plotpath 'dilationovertime.png'])
% saveas(h_forceovertime, [plotpath 'forceovertime.png'])

% saveas(h_vsdepth_forceanddil, [plotpath 'force_and_dilation_vs_time.png'])


% % partial
% findWindowPlots(dlength_partial,llength_partial,'partial dilation', 'magnitude', windowrange, wspacing, timeplot, ...
%     vfpdiff_colavg_partial, vfpdiff_colstd_partial, vfpdiff_rowavg_partial, vfpdiff_rowstd_partial, ...
%     vmag_colavg_partial,vmag_colstd_partial, vmag_rowavg_partial, vmag_rowstd_partial, ...
%     Fmag_colavg_partial, Fmag_colstd_partial, Fmag_rowavg_partial, Fmag_rowstd_partial, channelstring, grainstring,ynumgrids,znumgrids,Ltp)
%
%
% % compaction
% findWindowPlots(dlength_compact,llength_compact,'compaction', 'magnitude', windowrange, wspacing, timeplot, ...
%     vfpdiff_colavg_compact, vfpdiff_colstd_compact, vfpdiff_rowavg_compact, vfpdiff_rowstd_compact, ...
%     vmag_colavg_compact,vmag_colstd_compact, vmag_rowavg_compact, vmag_rowstd_compact, ...
%     Fmag_colavg_compact, Fmag_colstd_compact, Fmag_rowavg_compact, Fmag_rowstd_compact, channelstring, grainstring,ynumgrids,znumgrids,Ltp)
%
%
% % force components go in as absolute value so there is not missing data on
% % the log scale plots
% % --- x --- %
% %dilation
% findWindowPlots(dlength_dil,llength_dil,'dilation', 'X', windowrange, wspacing, timeplot, ...
%     vfpdiff_colavg, vfpdiff_colstd, vfpdiff_rowavg, vfpdiff_rowstd, ...
%     vx_colavg,vx_colstd, vx_rowavg, vx_rowstd, ...
%     abs(Fx_colavg), abs(Fx_colstd), abs(Fx_rowavg), abs(Fx_rowstd), channelstring, grainstring,ynumgrids,znumgrids,Ltp)
%
%
%
% % --- y --- %
% %dilation
% findWindowPlots(dlength_dil,llength_dil,'dilation', 'Y', windowrange, wspacing, timeplot, ...
%     vfpdiff_colavg, vfpdiff_colstd, vfpdiff_rowavg, vfpdiff_rowstd, ...
%     vy_colavg,vy_colstd, vy_rowavg, vy_rowstd, ...
%     abs(Fy_colavg), abs(Fy_colstd), abs(Fy_rowavg), abs(Fy_rowstd), channelstring, grainstring,ynumgrids,znumgrids,Ltp)
%
%
% % --- z --- %
% %dilation
% findWindowPlots(dlength_dil,llength_dil,'dilation', 'Z', windowrange, wspacing, timeplot, ...
%     vfpdiff_colavg, vfpdiff_colstd, vfpdiff_rowavg, vfpdiff_rowstd, ...
%     vz_colavg,vz_colstd, vz_rowavg, vz_rowstd, ...
%     abs(Fz_colavg), (Fz_colstd), abs(Fz_rowavg), abs(Fz_rowstd), channelstring, grainstring,ynumgrids,znumgrids,Ltp)






%% wave tracking at each depth

% the corresponding position in radial distance and depth is
% channel(yindex), depthplot(xindex)
Fmag_shocktrack = zeros(znumgrids,ynumgrids);
Fmag_shocktrack_onecross = zeros(znumgrids,ynumgrids);
Fmag_shocktrack_pos = zeros(znumgrids,ynumgrids);
Fmag_shocktrack_time = zeros(znumgrids,ynumgrids);
Fmag_shocktrack_time_onecross = zeros(znumgrids,ynumgrids);
Fmag_shocktrack_time_maxforce = zeros(znumgrids,ynumgrids);
vmag_shocktrack = zeros(znumgrids,ynumgrids);
vmag_shocktrack_time = zeros(znumgrids,ynumgrids);
Fmag_shocktrack_maxoverlap = zeros(znumgrids,ynumgrids);
Fmag_shocktrack_maxoverlap_time = zeros(znumgrids,ynumgrids);
Fmag_shocktrack_maxvelocity = zeros(znumgrids,ynumgrids);
Fmag_shocktrack_maxvelocity_time = zeros(znumgrids,ynumgrids);
Fmag_shocktrack_maxoverlap_first = zeros(znumgrids,ynumgrids);
Fmag_shocktrack_maxoverlap_first_time = zeros(znumgrids,ynumgrids);

loft_depth_track = zeros(znumgrids,ynumgrids);

rowcolorfull = turbo(znumgrids);

%threshold limits based on impact velocity
% if vel < 1% <= 0.1
%     shockstop = 10;%10;%20;
% else
%     shockstop = 30;%50;
% end


% set the bounds for computing the shock. do the first 50 m, index
%shockstop = 100; original

if SF
    shockstop = floor(SF*0.5/wspacing); %40;%20% 10;
else
    shockstop = floor(0.5/wspacing); %40;%20% 10;
end



shocktrigger = zeros(znumgrids,1);

% rowcolorfull = jet(length(1:13));
% figure
% axes('NextPlot','add')
% subplot(1,2,1)

segment_stop = round(3*ynumgrids/4);



for jj = 1:znumgrids%1:znumgridsFmag_logref_mean
%

    %h_testfig_forceovertime = figure;
    colplottemp = turbo(length(1:ynumgrids));

    for ii = 1:ynumgrids

%         if jj == 8 && ii > 34
%             [jj,ii]
%             keyboard
%         end

%         [Fmag_shocktrack_pos(jj,ii), max_force_ind] = max(store_Fmag{jj,ii});
%         Fmag_shocktrack_time_maxforce(jj,ii) = timeplot(max_force_ind);
        %mag_xy = (store_Fy{jj,ii}.^2+store_Fx{jj,ii}.^2).^0.5;
        %[Fmag_shocktrack(jj,ii) shockmaxdex] = max(mag_xy);

        % compute and store the stress magnitudes
%         store_psensor_mag_normal{jj,ii} = store_psensor_stress_xx{jj,ii}.*store_psensor_stress_xx{jj,ii} + ...
%                 store_psensor_stress_yy{jj,ii}.*store_psensor_stress_yy{jj,ii} + ...
%                 store_psensor_stress_zz{jj,ii}.*store_psensor_stress_zz{jj,ii};
%         store_psensor_mag_shear{jj,ii} = store_psensor_stress_xy{jj,ii}.*store_psensor_stress_xy{jj,ii} + ...
%                 store_psensor_stress_xz{jj,ii}.*store_psensor_stress_xz{jj,ii} + ...
%                 store_psensor_stress_yz{jj,ii}.*store_psensor_stress_yz{jj,ii};
%
%         store_mag_normal{jj,ii} = store_stress_xx{jj,ii}.*store_stress_xx{jj,ii} + ...
%                 store_stress_yy{jj,ii}.*store_stress_yy{jj,ii} + ...
%                 store_stress_zz{jj,ii}.*store_stress_zz{jj,ii};
%         store_mag_shear{jj,ii} = store_stress_xy{jj,ii}.*store_stress_xy{jj,ii} + ...
%                 store_stress_xz{jj,ii}.*store_stress_xz{jj,ii} + ...
%                 store_stress_yz{jj,ii}.*store_stress_yz{jj,ii};

        Fmag_p = store_Fmag{jj,ii};
        pdiffp = zeros(1,Ltp);
        pdiffp(2:Ltp) = (Fmag_p(2:Ltp))./Fmag_p(1:(Ltp-1));
        %[~, shockmaxdex] = max(pdiffp);

        %deriv = finite_diff_d1_fun(store_Fmag{jj,ii},wspacing);
        %[~, shockmaxdex] = max(deriv);


        % perform smoothing
        fw = 10;
        to_smooth = store_Fmag{jj,ii};
        y_smooth = smoothdata(to_smooth,'movmean',fw);
        deriv = finite_diff_d1_fun(y_smooth,wspacing);
        [~,derivmax] = max(deriv);
        [~,smoothmin] = min(y_smooth(1:derivmax));

        % new smoothing in parts
        % length of hi res data
%         L_hires = length(startFile:timestep1:endFile);
%         S_lowres = L_hires; % start index of low res data

%         y_smooth_hires = smoothdata(to_smooth(1:L_hires),'movmean',1); %8
%         x_smooth_hires =  timeplot(1:L_hires);
%
%         y_smooth_lores = smoothdata(to_smooth(S_lowres:Ltp),'movmean',1);
%         x_smooth_lores =  timeplot(S_lowres:Ltp);
%
%         % derivatives
%         dy_smooth_hires = finite_diff_d1_fun(y_smooth_hires,timestep1*simtimestep);
%         dy_smooth_lores = finite_diff_d1_fun(y_smooth_lores,timestep2*simtimestep);
%
%         % second derivatives
%         dy2_smooth_hires = finite_diff_d1_fun(dy_smooth_hires,timestep1*simtimestep);
%         dy2_smooth_lores = finite_diff_d1_fun(dy_smooth_lores,timestep2*simtimestep);
%
%
%         dy_combo = [dy_smooth_hires dy_smooth_lores(2:end)];
%         dy2_combo = [dy2_smooth_hires dy2_smooth_lores(2:end)];




        % find the first order of magnitude crossing using the second derivative data
        % if earlier in the channel, algorithm will select first point
        % everytime so just use first point for floor. if later in the
        % channel, use an average
%         if ii < 5
%             %floordy2dt2 = abs(dy2_combo(1));
%             %floordy2dt2 = store_Fmag{jj,ii}(1);
%             floordy2dt2 = mean(store_Fmag{jj,ii}(1:5));
%         else
%             %floordy2dt2 = mean(abs(dy2_combo(1:5)));
%             %floordy2dt2 = mean(store_Fmag{jj,ii}(1:5));
%             floordy2dt2 = mean(store_Fmag{jj,ii}(1:10));
%         end

        % threshold based on impact velocity
        if vel < 1 %0.01
            thresh = 0.15; %0.25;
            %thresh = 1;
        else
            thresh = 0.4;
        end
        thresh = 1;
        % minimum force in channel
%         [minforce,~] = min(y_smooth_hires);


        % find the minimum between 0 and the max force time, then find the
        % threshold crossing beyond that point

        [~, max_force_ind] = max(store_Fmag{jj,ii});
        %floor_force = Fmag_rowavg_init(1,jj);
        floor_force = 1E-8;

        mincheck = store_Fmag{jj,ii};
        mincheck(max_force_ind:Ltp) = 10E6; % set to some very large value
        [~,earlymin] = min(mincheck);

        threshcheck = store_Fmag{jj,ii};
        %threshcheck(1:earlymin) = 0; % set the earliest points to 0 so they aren't seen as the crossing

        %[Fmag_shocktrack_pos(jj,ii), max_force_ind] = max(store_Fmag{jj,ii});
%         [Fmag_shocktrack_pos(jj,ii), max_force_ind] = max(threshcheck);
%         Fmag_shocktrack_time_maxforce(jj,ii) = timeplot(max_force_ind);

%         neg1cross = find(abs(dy2_combo) > (1E-1), 1, 'first');
        %neg1cross = find(abs(dy2_combo) > ((10^(1.5))*floordy2dt2), 1, 'first');
        %neg1cross = find(store_Fmag{jj,ii} > ((10^(thresh))*floor_force), 1, 'first');

        cutoff_force = 10^(thresh)*floor_force;
        % check if this sensor has a higher starting force than the floor
        % if it doesn't, start doubling the threshold until it does
        numthresh = 1;
        while cutoff_force < store_Fmag{jj,ii}(1)
            cutoff_force = 10^(2*thresh)*store_Fmag{jj,ii}(1);
%             cutoff_force = 10^(numthresh*2*thresh)*floor_force;
            numthresh = numthresh+1;
        end

%         if jj == 8 && ii > 34
%             [jj,ii]
%             keyboard
%         end

        neg1cross = find(threshcheck > cutoff_force, 1, 'first');
        maxcheck = store_Fmag{jj,ii};
        maxcheck(max_force_ind:end) = 1000; % set high, so not slected
        %neg1cross = find(maxcheck < ((10^(thresh))*floor_force), 1, 'last');

        realmin = 1;
        if isempty(neg1cross)
            neg1cross = 1;
            %neg1cross = earlymin;
            %keyboard
            % update the shock trigger to be the previous value, if it's
            % the first time to update
            if shocktrigger(jj,1) == 0 && ii > 10
                shocktrigger(jj,1) = ii - 1;
            end
        else
            % step backward from this negative one crossing to find the "last
            % minimum" before the increase
            realmin = 1;
            realctr = 1;
            realctrold = neg1cross;
            minctrflag = 0;
            numcheck = 0;

            % check at most 20 steps backward
            while minctrflag < 1 && numcheck < neg1cross && realctr > 0 && realctrold > 1 %(ii - 5)
                realctr = realctrold - 1;
                [realctr, realctrold, jj, ii];
                 if (realctr < 1) || (realctrold < 1)
                     minctrflag = minctrflag + 1;
                 end

                newdiff = store_Fmag{jj,ii}(realctrold) - store_Fmag{jj,ii}(realctr);

                % if there we found the minimum OR the point when we
                % reached the noise floor, log the point

                if vel <= 0.5% 1
                    numdiv = 10;
                else
                    numdiv = 0.5;
                end
                %numdiv = 10;

                if newdiff < 0 || ...
                        newdiff < floor_force/numdiv
                    if  store_Fmag{jj,ii}(realctrold) < cutoff_force %floor_force*10^thresh
%                   if  store_Fmag{jj,ii}(realctrold) < minforce*10^1 % if on the same order of magnitude as the min
                        realmin = realctrold; % this is the index of first arrival
                        minctrflag = minctrflag + 1;
                    else
                        realctrold = realctrold - 1;
                        numcheck = numcheck + 1;
                    end
                else
                    % update the counters and termporary variables for the next
                    % loop
                    realctrold = realctrold - 1;
                    numcheck = numcheck + 1;
                end

            end

        end


%         if jj == 8 && ii > 34
%             [jj,ii]
%             keyboard
%         end

%         % overide if needed
%         [a_up,p_up,m_up] = force_pos_fun(jj,ii,comptype,vel);
%         if a_up
%             %Fmag_shocktrack(jj,ii) = store_Fmag{jj,ii}(a_up);
%             Fmag_shocktrack_time_onecross(jj,ii) = timeplot(a_up);
%             realmin = a_up;
%         end

        % log the time at index of first arrival
%         if realmin == 1
%             Fmag_shocktrack_time_onecross(jj,ii) = NaN;
%         else
%             Fmag_shocktrack_time_onecross(jj,ii) = timeplot(realmin);
%         end
        Fmag_shocktrack_time_onecross(jj,ii) = timeplot(neg1cross);
        %Fmag_shocktrack_time_onecross(jj,ii) = timeplot(realmin);
        % max force storing - check for within some time of first arrival,
        % based on depth
        maxcheck = store_Fmag{jj,ii};
        % reset values outside range (below arrival, above
        % arrival+threshold) to zero
        % timeplot increment is 0.001 with hi res data
        % looking 0.05 s ahead = 50 time steps
        lookahead = 15;



        [Fmag_shocktrack_pos(jj,ii), max_force_ind] = max(threshcheck);
        Fmag_shocktrack_time_maxforce(jj,ii) = timeplot(max_force_ind);

        maxcheck(1:(realmin-1)) = 0;
        maxcheck((realmin+lookahead):end) = 0;


        % this is the "virtual peak force", use it to make sure the peak
        % force identified is on the order of the peak force in the cell
        [virtual_peak_force, peak_force_ind] = max(maxcheck);




        % find the 0 crossing after the first maximum from derivative data
        % - this is the peak force
        % use negative one crossing instead of finding the first peak
%         [temppeaks,tempfind] = findpeaks(dy_combo);
        % temppeaks is the force value of the peaks
        % tempfind is the index of the peak forces

        % want to find the peaks after the negative one crossing
        % get a list of them, then check through them iteratively until the
        % identified peak is within an order of magnitude of the peak force
        % ind value

        % tempfind_plot
%         tempfind_plot = find(tempfind >= neg1cross, 1, 'first');
        % find the first peak that is greater than the derivative crossing

        %keyboard

        tempfind = [];
        if isempty(tempfind)
            tempfind = 5;
        end
        tempfind_plot = [];
        if isempty(tempfind_plot)
            tempfind_plot = 1;
        end

        % set threshold based on impacct velocity
        if vel < 10^-1
            derivthresh = -0.5;%-0.14;
        elseif vel < 10^0
            derivthresh = -0.5;
        elseif vel == 1
            derivthresh = -0.75;
        else
            derivthresh = -1;
        end

        pff = 0; % peak found flag
        p_ctr = 1; % peak counter to exit while loop if stuck

%         while (pff < 1) && (tempfind_plot < length(tempfind))

            % zero out the values up to the temporary peak  we are at
%             %dy_combo(tempfind(1:tempfind_plot)) = 0;
%             dy_combo_temp = dy_combo;
%             %dy_combo_temp(1:tempfind(tempfind_plot)) = 0;
%             dy_combo_temp(1:realmin) = 0;
%             %smoothdata(to_smooth(1:L_hires),'movmean',1);
%
%             % this finds the first time the derivative crosses 0 after the
%             % peak selected
%             zerocross = find(dy_combo_temp < 0, 1, 'first')-1;
%
            % if this peak is not within a threshold value of the max peak,
            % reset tempfind_plot to the next index and try again
%             if store_Fmag{jj,ii}(zerocross) > (10^derivthresh)*virtual_peak_force
%                 pff = pff + 1; % exit the loop
%             else
%                 tempfind_plot = tempfind_plot + 1;
%                 p_ctr = p_ctr + 1;
%             end
%
%         end % end zero cross loop

%         % if it didn't get trigger once in the loop, generate an empty
%         % value
%         if pff == 0
%             dy_combo_temp = dy_combo;
%             dy_combo_temp(1:tempfind(tempfind_plot)) = 0;
%             zerocross = find(dy_combo_temp < 0, 1, 'first')-1;
%         end


%         if jj == 8 && ii > 34
%             [jj,ii]
%             keyboard
%         end


        zerocross = [];
        if isempty(zerocross)
            zerocross = Ltp - 10;
        end
        if zerocross <= 0
            zerocross = 1;
        end


        % compute the order of magnitude of each force value above the
        % minimum force value for this cell
        [shock_min_F, shockmindex] = min(store_Fmag{jj,ii});
        oom_above_min = log10(store_Fmag{jj,ii}./shock_min_F);
        % find the first instance of crossing 1
        [~,min_oom_index] = min(oom_above_min);
        %oom_above_min(1:min_oom_index) = 0;
        oom_above_min(1:smoothmin) = 0;
        above1 = oom_above_min - 1;
        onecrossing = find(above1 > 0, 1, 'first');



        % zero out values before the minimum index in case the force in the
        % cell started higher than the peak force seen due to the shock
        %shockmaxdex = onecrossing;
        % find the index when t = 0.2 s, then find the minimum point in
        % that range
        temp_tp = round(timeplot*100);
        p2_min = find(temp_tp == 20);
        y_eval = y_smooth;
        y_eval(1:smoothmin) = 0;
%         y_eval(1:min_oom_index) = 0;
        %y_eval(p2_min:Ltp) = 0;
        [~,shockmaxdex] = max(y_eval);

        %keyboard
        % store one crossing for first arrival
        if isempty(onecrossing) onecrossing = Ltp; end

        if vel < 1 % for low velocities
            if onecrossing < 4
                one_ind = Ltp;
            else
                one_ind = onecrossing - 3;
            end
        else
            one_ind = onecrossing;
        end



        % reset to be from the derivative zero crossing method
         shockmaxdex = zerocross;
   %     shockmaxdex = peak_force_ind;

        if isempty(shockmaxdex)
            %keyboard
            shockmaxdex = Ltp;
        end


        % find an exponential fit between the min and the zerocrossing
%         if zerocross <= realmin
%             keyboard
%         end

        if ii <= shockstop
            if zerocross <= realmin
                %keyboard
                % look ahead 10 time steps if the zerocross didnt' work
                zerocross = realmin + 5;
            end

            % setup exponential function and input data
            % add a condition to shorten the exp_range fit to only the
            % first few points after the real min if we are at early radial
            % positions
            first_fit_range = realmin:zerocross;
            if ii < 10
                first_fit_range = realmin:(realmin+2);
            end

%             exp_fit_time = timeplot(first_fit_range).';
%             exp_fit_y = store_Fmag{jj,ii}(first_fit_range);
%             myexpfun = @(x,xdata) x(1)*exp(x(2)*xdata);
%             % calculate initial guess for exp model based off endpoints
%             % k growth rate, A amplitude
%             k = (1/(exp_fit_time(1) - exp_fit_time(end)))*log(exp_fit_y(1)/exp_fit_y(end));
%             A = exp_fit_y(1)*exp(-1*k*exp_fit_time(1));
%             %p0_exp = [floor_force,1000];
%             p0_exp = [A,k];
%             p_exp = lsqcurvefit(myexpfun,p0_exp,exp_fit_time,exp_fit_y);
%             A_rise = p_exp(1);
%             k_rise = p_exp(2);
%             %[jj,ii]
%             % now extend the fit to reach the max value
%             % find the time to use (time exponential would reach max)
%             tmax = log((store_Fmag{jj,ii}(max_force_ind)/A))/k;
%             % add a little padding
%
%             t_exp_ind = find(timeplot < tmax,1,'last')+5;
%             if isempty(t_exp_ind)
%                 t_exp_ind = Ltp;
%             end
%             if t_exp_ind > Ltp
%                 t_exp_ind = Ltp;
%             end
%             exp_range = realmin:t_exp_ind;
%
%             new_exp_time = timeplot(exp_range).';
%
%             exp_fit_force = myexpfun(p_exp,new_exp_time);
% %             exp_fit_force = myexpfun(p_exp,exp_fit_time);
%
%             % find the percent difference between the exponential fit
%             % points and the data, the last point within #% of fit is the
%             % ending exponential point
% %             if evalval <= 1 && ii > 15
% %                 pdiff_allow = -0.9;
% %             else
% %                 pdiff_allow = -0.5;
% %             end
%             pdiff_allow = -0.5; % should be the same for all tests, otherwise tracking wrong point
%             exp_data_comp = store_Fmag{jj,ii}(exp_range);
%             % percent diff between data and fit
%             exp_pdiff = (exp_data_comp - exp_fit_force)./exp_fit_force;
%             last_exp_ind = find(exp_pdiff > pdiff_allow,1,'last');
%             actual_exp_time_ind = find(timeplot == new_exp_time(last_exp_ind));
%
%
%             % now fit the ending data
%             decay_range = (t_exp_ind-5):Ltp;
%             % reset to the decay range to only capture the first "ring" of
%             % the ringing portion for high speed piston impacts
%             if evalval >= 2
%                 if t_exp_ind+30 < Ltp
%                     decay_range = (t_exp_ind-5):(t_exp_ind+30);
%                 end
%             end
%             if isempty(decay_range)
%                 keyboard
%                 % use half the data if the zero crossing had an issue
%                 decay_range = round(Ltp/2):Ltp;
%             end
%             %keyboard
%             decay_time = timeplot(decay_range).';


%             [time_uniform_decay,gen_uniform_decay,A_decay,k_decay] =  ...
%                 least_fit_decay_fun_overlaps(decay_time,store_Fmag{jj,ii}(decay_range),evalval,jj);
%
%
%             % find the intersection of the exponential functions
%             t_intersect = log(A_decay/A_rise)/(-k_decay + k_rise);
%             F_intersect = myexpfun(p_exp,t_intersect);



        end


%         if jj == 8 && ii > 34
%             [jj,ii]
%             keyboard
%         end

%% old max stress method of peak finding
            [maxstress,max_stresses_inds] = findpeaks(-1*store_stress_yy{jj,ii});

            %[~,max_stress_ind] = min(store_stress_yy{jj,ii})
            % assume the peak is the first found, but check if it is at
            % least 1% greater than starting value. if not, go to next
            max_stress_ind = 1;
            max_stress_count = 1;
            max_stress_tot = length(max_stresses_inds);
            max_stress_foundflag = 0;
            init_yy_stress = -1*store_stress_yy{jj,ii}(1);
            while max_stress_foundflag < 1 && max_stress_ind < max_stress_tot
                check_stress = 100*(-1*store_stress_yy{jj,ii}(max_stresses_inds(max_stress_ind)) - ...
                    init_yy_stress)/init_yy_stress;
                if abs(check_stress) > 0.015 %0.01
                    max_stress_foundflag = max_stress_foundflag+1;
                else
                    max_stress_ind = max_stress_ind + 1;
                end
            end

%         if jj == 8 && ii > 34
%             [jj,ii]
%             keyboard
%         end
%
%             if ~isempty(max_stresses_inds)
%                 Fmag_shocktrack_time(jj,ii) = timeplot(max_stresses_inds(max_stress_ind));
%                 Fmag_shocktrack(jj,ii) = store_Fmag{jj,ii}(max_stresses_inds(max_stress_ind));
%             else
%                 Fmag_shocktrack_time(jj,ii) = nan;
%                 Fmag_shocktrack(jj,ii) = nan;
%             end






    %% new delta method of peak finding
%             if jj > 9
%                 keyboard
%             end
            [maxdelta,max_deltas_inds] = findpeaks(store_delta{jj,ii});
            [maxdelta_max,max_deltas_inds_max] = max(store_delta{jj,ii});

            %[~,max_stress_ind] = min(store_stress_yy{jj,ii})
            % assume the peak is the first found, but check if it is at
            % least 1% greater than starting value. if not, go to next
            max_delta_ind = 1;
            max_delta_count = 1;
            max_delta_tot = length(max_deltas_inds);
            max_delta_foundflag = 0;
            init_delta = store_delta{jj,ii}(1);
            while max_delta_foundflag < 1 && max_delta_ind < max_delta_tot
                check_delta = 100*(store_delta{jj,ii}(max_deltas_inds(max_delta_ind)) - ...
                    init_delta)/init_delta;
                if abs(check_delta) > 0.015 %0.01
                    max_delta_foundflag = max_delta_foundflag+1;
                else
                    max_delta_ind = max_delta_ind + 1;
                end
            end

            Fmag_shocktrack_maxoverlap(jj,ii) = maxdelta_max;
            Fmag_shocktrack_maxoverlap_time(jj,ii) = timeplot(max_deltas_inds_max);
%
%
%
%             if ~isempty(max_deltas_inds)
%                 Fmag_shocktrack_time(jj,ii) = timeplot(max_deltas_inds(max_delta_ind));
%                 Fmag_shocktrack(jj,ii) = store_Fmag{jj,ii}(max_deltas_inds(max_delta_ind));
%             else
%                 Fmag_shocktrack_time(jj,ii) = nan;
%                 Fmag_shocktrack(jj,ii) = nan;
%             end
%
%
%
%
%     [maxyforce,max_yforce_inds] = max(store_Fy{jj,ii});



 %%
             [maxyforce,max_yforce_inds] = findpeaks(abs(store_Fy{jj,ii}));

            %[~,max_stress_ind] = min(store_stress_yy{jj,ii})
            % assume the peak is the first found, but check if it is at
            % least 1% greater than starting value. if not, go to next
            max_yforce_ind = 1;
            max_yforce_count = 1;
            max_yforce_tot = length(max_yforce_inds);
            max_yforce_foundflag = 0;

            % want to take an average of the first two points to get the
            % initial force to average out any noise. however,
            % when the cell is too close to the piston cite there is not enough time for the
            % force to evolve before the peak arrives, so for the early
            % cells, just use the first value
            % tuned to use less cells at top speed



            if vel < 10
                limitercell = 6;
            else
                limitercell = 7;
            end


            if ii < limitercell
                init_yforce = mean(abs(store_Fy{jj,ii}(1)));
            else
                init_yforce = mean(abs(store_Fy{jj,ii}(1:2)));
            end

            % 15 heuristically chosen
            if ii > 15 && vel < 10
                init_yforce = mean(abs(store_Fy{jj,ii}(1:4)));
            end

            if ii > 35 && vel == 10
                init_yforce = mean(abs(store_Fy{jj,ii}(1:8)));
            end

            if ii > 45 && vel == 10
                init_yforce = mean(abs(store_Fy{jj,ii}(1:15)));
            end

            if ii > 65 && vel == 10
                init_yforce = mean(abs(store_Fy{jj,ii}(1:25)));
            end

            thresh_val = 10;
            if  vel == 10
                thresh_val = 100;
            end

            if r_p == 5E-3 || r_p == 12.5E-3 || r_p == 0.5E-3
                thresh_val = 100000;
            end

            %keyboard


            while max_yforce_foundflag < 1 && max_yforce_ind < max_yforce_tot
                check_yforce = (abs(store_Fy{jj,ii}(max_yforce_inds(max_yforce_ind))) - ...
                    init_yforce)/init_yforce;
                if abs(check_yforce) > thresh_val %0.01
                    max_yforce_foundflag = max_yforce_foundflag+1;
                else
                    max_yforce_ind = max_yforce_ind + 1;
                end
            end



            if ~isempty(max_yforce_inds)
                Fmag_shocktrack_time(jj,ii) = timeplot(max_yforce_inds(max_yforce_ind));
                Fmag_shocktrack(jj,ii) = store_Fy{jj,ii}(max_yforce_inds(max_yforce_ind));

                Fmag_shocktrack_maxoverlap_first(jj,ii) = store_delta{jj,ii}(max_yforce_ind);
                Fmag_shocktrack_maxoverlap_first_time(jj,ii) = timeplot(max_yforce_ind);

                % now check to see if we are still on exponential rise, if
                % not, step back until we are
                if ii > range_exp && ii < range_prop

                    pdiff_last = 0.001;
                    c_ind = max_yforce_inds(max_yforce_ind);


                                        % uncomment below for the height sweep case
%                     if comptype == '2m20cm_Pwave_loose_velocitysweep'
%                         comptype_ind_flag = 1;
%                     else
%                         comptype_ind_flag = 0;
%                     end
%
%                     if (comptype_ind_flag == 1) && vel > 0.6
%                         limiter_percent = 2;
%                     else
%                         limiter_percent = 1;
%                     end
                    limiter_percent = 1;

                    % this was used for finding the 'correct' peak when
                    % sound waves were involved, but in paper 2 we don't do
                    % that, so let's turn off since it is presenting issues
%                     while pdiff_last < limiter_percent %1
%
%                         pdiff_last = (abs(store_Fy{jj,ii}(c_ind)) - abs(store_Fy{jj,ii}(c_ind-1)))/abs(store_Fy{jj,ii}(c_ind-1));
%                         c_ind = c_ind - 1;
%                         if c_ind < 2
%                             c_ind = 2;
%                             keyboard
%                         end
%                     end

                    %max_yforce_ind = c_ind;
                    c_ind = c_ind + 1; % correct to the second to last step
                    Fmag_shocktrack_time(jj,ii) = timeplot(c_ind);
                    Fmag_shocktrack(jj,ii) = store_Fy{jj,ii}(c_ind);

                    Fmag_shocktrack_maxoverlap_first(jj,ii) = store_delta{jj,ii}(c_ind);
                    Fmag_shocktrack_maxoverlap_first_time(jj,ii) = timeplot(c_ind);

                end

            else
                Fmag_shocktrack_time(jj,ii) = nan;
                Fmag_shocktrack(jj,ii) = nan;

                Fmag_shocktrack_maxoverlap_first(jj,ii) = nan;
                Fmag_shocktrack_maxoverlap_first_time(jj,ii) = nan;
            end

    %%
%             [maxyforce,max_yforce_inds] = findpeaks(store_Fy{jj,ii});
%
%             %[~,max_stress_ind] = min(store_stress_yy{jj,ii})
%             % assume the peak is the first found, but check if it is at
%             % least 1% greater than starting value. if not, go to next
%             max_yforce_ind = 1;
%             max_yforce_count = 1;
%             max_yforce_tot = length(max_yforce_inds);
%             max_yforce_foundflag = 0;
%             init_yforce = store_Fy{jj,ii}(1);
%             while max_yforce_foundflag < 1 && max_yforce_ind < max_yforce_tot
%                 check_yforce = 100*(store_Fy{jj,ii}(max_yforce_inds(max_yforce_ind)) - ...
%                     init_yforce)/init_yforce;
%                 if abs(check_yforce) > 0.01
%                     max_yforce_foundflag = max_yforce_foundflag+1;
%                 else
%                     max_yforce_ind = max_yforce_ind + 1;
%                 end
%             end
%
%             Fmag_shocktrack_time(jj,ii) = timeplot(max_yforce_inds(max_yforce_ind));
%             Fmag_shocktrack(jj,ii) = store_Fmag{jj,ii}(max_yforce_inds(max_yforce_ind));
%             Fmag_shocktrack_time(jj,ii) = timeplot(max_yforce_inds);
%             Fmag_shocktrack(jj,ii) = store_Fmag{jj,ii}(max_yforce_inds);
%             if t_intersect < 0


%             if t_intersect < 0
%                 % go plot stuff to see what's wrong
%                 keyboard
%             end


%         else
%             Fmag_shocktrack(jj,ii) = store_Fmag{jj,ii}(zerocross);
%             Fmag_shocktrack_time(jj,ii) = timeplot(zerocross);
%         end
%         Fmag_shocktrack(jj,ii) = store_Fmag{jj,ii}(actual_exp_time_ind);
        %deriv = finite_diff_d1_fun(store_Fmag{jj,ii},wspacing);
%         Fmag_shocktrack_time(jj,ii) = timeplot(shockmaxdex);
        %Fmag_shocktrack_time(jj,ii) = timeplot(shockmaxdex);
        [vmag_shocktrack(jj,ii), shockmaxdex_v] = max(store_Vmag{jj,ii});
        vmag_shocktrack_time(jj,ii) = timeplot(shockmaxdex_v);


        % overwrite peak and arrival if manual update specified
%         [a_up,p_up,m_up] = force_pos_fun(jj,ii,comptype,evalval);
        %[a_up,p_up,m_up] = force_pos_fun_1ms_tracking_wrongpeak_02142022(jj,ii,comptype,evalval);
        % if either index update is specified, update stored time for this
        % cell
%         if a_up
%             %Fmag_shocktrack(jj,ii) = store_Fmag{jj,ii}(a_up);
%             Fmag_shocktrack_time_onecross(jj,ii) = timeplot(a_up);
%             realmin = a_up;
%         end
        if 0%p_up
            Fmag_shocktrack(jj,ii) = store_Fmag{jj,ii}(p_up);
            Fmag_shocktrack_time(jj,ii) = timeplot(p_up);
            peak_force_ind = p_up;
        end
        if 0%m_up
            Fmag_shocktrack_pos(jj,ii) = store_Fmag{jj,ii}(m_up);
            Fmag_shocktrack_time_maxforce(jj,ii) = timeplot(m_up);
            max_force_ind = m_up;
        end


        if jj == sz && ii == sr
            sensor_force = store_Fy{jj,ii};
            sensor_force_arrival = realmin;
            sensor_force_peak = Fmag_shocktrack(jj,ii);
            sensor_force_peak_time = Fmag_shocktrack_time(jj,ii);

            sensor_stressyy = store_stress_yy{jj,ii};
%             sensor_stress_peak_time = timeplot(max_stresses_inds(max_stress_ind));
            sensor_stress_peak_time = nan;
            %sensor_stress_peak_val = store_stress_yy{jj,ii}(max_stresses_inds(max_stress_ind));
            sensor_stress_peak_val = nan;

            sensor_stress_yz = store_stress_yz{jj,ii};

            sensor_overlap = store_delta{jj,ii};
            sensor_overlap_peak = maxdelta_max;
            sensor_overlap_peak_time = timeplot(max_deltas_inds_max);

        end

        % store maximum particle velocity
        plot_Vmag = (store_VX{jj,ii}.*store_VX{jj,ii} + store_VY{jj,ii}.*store_VY{jj,ii} + store_VZ{jj,ii}.*store_VZ{jj,ii}).^(0.5);
        % store max vel
        [max_velmag,max_velmag_ind] = max(plot_Vmag);
        Fmag_shocktrack_maxvelocity(jj,ii) = max_velmag;
        Fmag_shocktrack_maxvelocity_time(jj,ii) = timeplot(max_velmag_ind);

%% lofted depth calculation
        % store 'lofted depth or not'
        normal_stress_mag = (store_stress_yy{jj,ii}.*store_stress_yy{jj,ii} + ...
        store_stress_xx{jj,ii}.*store_stress_xx{jj,ii} + ...
        store_stress_zz{jj,ii}.*store_stress_zz{jj,ii}).^(0.5);

        % find the minimum normal stress
        init_normal_stress = normal_stress_mag(1);
        min_normal_stress = min(normal_stress_mag);
        normal_stress_ratio = init_normal_stress/min_normal_stress;
        if normal_stress_ratio > 100
            loft_depth_track(jj,ii) = 1; % flag this cell
        end



        depthprint = 2;
        if  0 %jj == 5 %jj > 36 %jj > 10 %jj == 34 %jj == 17 || jj == 20 %jj == 30 || jj == 35 % || jj == 15
             %   || jj == 18
            if  0 % ii == 200 % <= shockstop %ii == 5 || ii == 15 || ii == 25 || ((ii < shockstop) && (ii > 60)) % || (ii==shockstop)


            %deriv = finite_diff_d1_fun(finite_diff_d1_fun(finite_diff_d1_fun(y_smooth,fw),fw),fw);
            %deriv = y_smooth; %smoothdata(y_smooth,'movmean',fw);
            figure
            subplot(1,1,1)
            semilogy(timeplot,store_Fmag{jj,ii},'k.','LineWidth',2,'DisplayName','Force Magnitude','LineWidth',2)
            %semilogy(timeplot,y_smooth,'k--','DisplayName','SingleSmooth','LineWidth',2)
            hold on
            grid on
            semilogy(timeplot(max_force_ind), store_Fmag{jj,ii}(max_force_ind),'rx','LineWidth',2,'DisplayName','Max Force','MarkerSize',10)
            %semilogy(x_smooth_hires,y_smooth_hires,'r--','DisplayName','HiRes','LineWidth',2)
            %semilogy(x_smooth_lores,y_smooth_lores,'b--','DisplayName','LoRes','LineWidth',2)
            semilogy(timeplot(zerocross),store_Fmag{jj,ii}(zerocross),'bo','DisplayName','ZeroCross','LineWidth',2)
            semilogy(timeplot(realmin),store_Fmag{jj,ii}(realmin),'ro','DisplayName','First Arrival','LineWidth',2,'MarkerSize',10)
            %semilogy(timeplot(neg1cross),store_Fmag{jj,ii}(neg1cross),'ro','DisplayName','First Arrival old','LineWidth',2)
%             semilogy(timeplot(peak_force_ind),store_Fmag{jj,ii}(peak_force_ind),'co','DisplayName','Peak Force','LineWidth',2,'MarkerSize',10)
            semilogy(timeplot,cutoff_force*ones(1,length(timeplot)),'r-','DisplayName','FloorForce*thresh')
           % semilogy(new_exp_time,exp_fit_force,'c--','DisplayName','Exp fit','LineWidth',1.5)
            %semilogy(time_uniform_decay,gen_uniform_decay,'g--','DisplayName','Decay fit','LineWidth',1.5)
            %semilogy(exp_fit_time,exp_fit_y,'bo','DisplayName','Points for exp fit')
            %semilogy(t_intersect,F_intersect,'r+','DisplayName','Intersection','MarkerSize',10,'LineWidth',1.5)
           % semilogy(timeplot(actual_exp_time_ind),store_Fmag{jj,ii}(actual_exp_time_ind),'co','DisplayName','Peak Force','LineWidth',2,'MarkerSize',10)


            title(['Force v - time, ' ', depth ind ' num2str(jj) ', radial ind ' num2str(ii)  ', v_i = ' num2str(vel)])
            legend('show','location','best')
            set(gca,'fontsize',14)
            xlabel('Time [s]')
            ylabel('Average Particle Force [N]')








            %figure - transverse stresses - off for now
%             subplot(1,3,3)
%             figure
%             semilogy(timeplot,abs(store_stress_xy{jj,ii}),'k.--','LineWidth',2,'DisplayName','xy','LineWidth',2)
%             %semilogy(timeplot,y_smooth,'k--','DisplayName','SingleSmooth','LineWidth',2)
%             hold on
%             grid on
%             semilogy(timeplot,abs(store_stress_xz{jj,ii}),'r.--','LineWidth',2,'DisplayName','xz','LineWidth',2)
%             semilogy(timeplot,abs(store_stress_yz{jj,ii}),'b.--','LineWidth',2,'DisplayName','yz','LineWidth',2)
%             title(['Shear stress, depth ind ' num2str(jj) ', radial ind ' num2str(ii)  ', v_i = ' num2str(vel)])
%             legend('show','location','best')
%             set(gca,'fontsize',14)
%             xlabel('Time [s]')
%             ylabel('Average Particle shear stress  [N/m/m]')
%
%             figure
%             plot(timeplot,store_vf_mean_post{jj,ii},'k.','LineWidth',2,'DisplayName','xy','LineWidth',2)
%             %semilogy(timeplot,y_smooth,'k--','DisplayName','SingleSmooth','LineWidth',2)
%             hold on
%             grid on
% %             semilogy(timeplot,store_stress_xz{jj,ii},'r.','LineWidth',2,'DisplayName','xz','LineWidth',2)
% %             semilogy(timeplot,store_stress_yz{jj,ii},'b.','LineWidth',2,'DisplayName','yz','LineWidth',2)
%             title(['\phi, depth ind ' num2str(jj) ', radial ind ' num2str(ii)  ', v_i = ' num2str(vel)])
%             legend('show','location','best')
%             set(gca,'fontsize',14)
%             xlabel('Time [s]')
%             ylabel('\phi')
%
%             figure
%             plot(timeplot,store_psensor_stress_yz{jj,ii},'k--.')
%             xlabel('time')
%             ylabel('yz shear')
%             title(['depth ind ' num2str(jj) ', radial ind ' num2str(ii) ])
%             set(gca,'fontsize',14)
%
%             figure
%             plot(timeplot,store_psensor_stress_yy{jj,ii},'k--.')
%             xlabel('time')
%             ylabel('yy normal')
%             title(['depth ind ' num2str(jj) ', radial ind ' num2str(ii) ])
%             set(gca,'fontsize',14)

%             figure
%             semilogy(timeplot,store_psensor_delta{jj,ii},'k--.')
%             xlabel('time')
%             ylabel('\delta [m]')
%             title(['depth ind ' num2str(jj) ', radial ind ' num2str(ii) ])
%             set(gca,'fontsize',14)


%             figure
% %             stress_mag = store_stress_xy{jj,ii}.*store_stress_xy{jj,ii} + ...
% %                 store_stress_xz{jj,ii}.*store_stress_xz{jj,ii} + ...
% %                 store_stress_yz{jj,ii}.*store_stress_yz{jj,ii};
% %            semilogy(timeplot,store_mag_shear{jj,ii},'k--.','DisplayName','Sensor Shear Mag')
%             grid on
%             hold on
%             semilogy(timeplot,abs(store_stress_yy{jj,ii}),'r--.','DisplayName','Sensor Normal')
%             grid on
%             hold on
%             %semilogy(timeplot,store_psensor_mag_shear{jj,ii},'b--.','DisplayName','Particle Shear Mag')
%             semilogy(timeplot,abs(store_psensor_stress_yy{jj,ii}),'c--.','DisplayName','Particle Normal')
%             %semilogy(timeplot,store_psensor_delta{jj,ii},'g--.','DisplayName','Delta')
%             xlabel('time')
%             ylabel('shear stress magnitude [N/m/m]')
%             title(['depth ind ' num2str(jj) ', radial ind ' num2str(ii) ])
%             set(gca,'fontsize',14)
%             legend('show','location','best')
% %             xlim([0, 0.06])

%% normal stress plot

            normal_stress_mag = (store_stress_yy{jj,ii}.*store_stress_yy{jj,ii} + ...
            store_stress_xx{jj,ii}.*store_stress_xx{jj,ii} + ...
            store_stress_zz{jj,ii}.*store_stress_zz{jj,ii}).^(0.5);

            figure
            plot(timeplot,abs(store_stress_yy{jj,ii}),'r--.','DisplayName','YY')
            grid on
            hold on
            plot(timeplot,abs(store_stress_xx{jj,ii}),'b--.','DisplayName','XX')
            plot(timeplot,abs(store_stress_zz{jj,ii}),'k--.','DisplayName','ZZ')
            plot(timeplot,normal_stress_mag,'c--.','DisplayName','magnitude')
            xlabel('time s')
            legend('show','location','best')
            ylabel('|normal stress|)')
            xlim([timeplot(90), timeplot(230)])
            title(['depth ' num2str(depthplot(jj)) ' cm'])

%% multi stress plot
            figure
            subplot(1,2,1)
            plot(timeplot,(store_stress_yz{jj,ii}),'r--.','DisplayName','YZ')
            grid on
            hold on
            plot(timeplot,(store_stress_xy{jj,ii}),'b--.','DisplayName','XY')
            plot(timeplot,(store_stress_xz{jj,ii}),'k--.','DisplayName','YZ')
            xlabel('time s')
            legend('show','location','best')
            ylabel('shear stress')
            xlim([timeplot(90), timeplot(230)])
            title(['depth ' num2str(depthplot(jj)) ' cm'])

            %figure
            subplot(1,2,2)
            plot(timeplot,(store_stress_yy{jj,ii}),'r--.','DisplayName','YY')
            grid on
            hold on
            plot(timeplot,(store_stress_xx{jj,ii}),'b--.','DisplayName','XX')
            plot(timeplot,(store_stress_zz{jj,ii}),'k--.','DisplayName','zZ')
            xlabel('time s')
            ylabel('normal stress')
            legend('show','location','best')
            xlim([timeplot(90), timeplot(230)])
            title(['depth ' num2str(depthplot(jj)) ' cm'])

            figure
            %subplot(1,2,1)
            semilogy(timeplot,store_delta{jj,ii},'k.--','LineWidth',2,'DisplayName','Average','LineWidth',2)
            %semilogy(timeplot,y_smooth,'k--','DisplayName','SingleSmooth','LineWidth',2)
            hold on
            grid on
            %semilogy(timeplot,store_psensor_delta{jj,ii},'r--.','DisplayName','Sensor')
            semilogy(timeplot(max_deltas_inds(max_delta_ind)),store_delta{jj,ii}(max_deltas_inds(max_delta_ind)),'go','LineWidth',2,'DisplayName','Max','LineWidth',2)

            title(['Overlap v - time, ' ', depth ind ' num2str(jj) ', radial ind ' num2str(ii)  ', v_i = ' num2str(vel)])
            legend('show','location','best')
            set(gca,'fontsize',14)
            xlabel('Time [s]')
            ylabel('\delta [m]')

 %% return

                        %figure
%             %subplot(1,2,2)
%             figure
%             xx_plot = store_stress_xx{jj,ii};
%             xx_plot(xx_plot>= 0) = nan;
%             yy_plot = store_stress_yy{jj,ii};
%             yy_plot(yy_plot>= 0) = nan;
%             zz_plot = store_stress_zz{jj,ii};
%             zz_plot(zz_plot>= 0) = nan;
%             semilogy(timeplot,xx_plot,'k.--','LineWidth',2,'DisplayName','xx','LineWidth',2)
%             %semilogy(timeplot,y_smooth,'k--','DisplayName','SingleSmooth','LineWidth',2)
%             hold on
%             grid on
%             semilogy(timeplot,yy_plot,'r.--','LineWidth',2,'DisplayName','yy','LineWidth',2)
%             semilogy(timeplot,zz_plot,'b.--','LineWidth',2,'DisplayName','zz','LineWidth',2)
%             semilogy(timeplot(max_stresses_inds(max_stress_ind)),store_stress_yy{jj,ii}(max_stresses_inds(max_stress_ind)),'go','LineWidth',2,'DisplayName','Peak yy','LineWidth',2)
%             title(['Stress - time, ' ', depth ind ' num2str(jj) ', radial ind ' num2str(ii)  ', v_i = ' num2str(vel)])
%             legend('show','location','best')
%             set(gca,'fontsize',14)
%             xlabel('Time [s]')
%             ylabel('Average Particle normal stress  [N/m/m]')


            figure
            subplot(1,1,1)
            semilogy(timeplot,abs(store_Fy{jj,ii}),'k--.','LineWidth',1.5,'DisplayName','Cell Censor','LineWidth',2)
            grid on
            hold on
            semilogy(timeplot(max_yforce_inds(max_yforce_ind)),abs(store_Fy{jj,ii}(max_yforce_inds(max_yforce_ind))),'ro','LineWidth',2,'DisplayName','Force Y','LineWidth',2)
            %semilogy(timeplot,abs(store_psensor_Fy{jj,ii}),'r--.','LineWidth',1.5,'DisplayName','Particle Sensor','LineWidth',2)
            legend('show','location','best')
            set(gca,'fontsize',14)
            xlabel('Time [s]')
            ylabel('Average Particle Force [N]')
            xlim([0,0.05])
            title(['Fy - time, ' ', depth ind ' num2str(jj) ', radial ind ' num2str(ii)  ', v_i = ' num2str(vel)])

            %se

%             figure
%             semilogy(timeplot,store_psensor_mag_normal{jj,ii},'k--.')
%             xlabel('time')
%             ylabel('normal stress magnitude [N/m/m]')
%             title(['depth ind ' num2str(jj) ', radial ind ' num2str(ii) ])
%             set(gca,'fontsize',14)
        keyboard



%             subplot(3,1,3)
%             semilogy(x_smooth_hires,abs(dy2_smooth_hires),'r-','DisplayName','HiRes')
%             hold on
%             grid on
%             semilogy(x_smooth_lores,abs(dy2_smooth_lores),'b-','DisplayName','LoRes')
%             semilogy(timeplot(neg1cross),abs(dy2_combo(neg1cross)),'ro','LineWidth',2,'DisplayName','10^{-1} Crossing')
%             title('Second Derivatives')
%             legend('show','location','best')

        %keyboard
%
            %pdiff plot
%             Fmag_p = store_Fmag{jj,ii};
%             pdiffp = zeros(1,Ltp);
%             pdiffp(2:Ltp) = (Fmag_p(2:Ltp))./Fmag_p(1:(Ltp-1));
%
%             figure
%             semilogy(timeplot,pdiffp,'k--.')
%             hold on
%             title('percent difference from prior point')


%             figure
%             % percent difference above min
%             [shock_min_F, shockmindex] = min(store_Fmag{jj,ii});
%             oom_above_min = log10(store_Fmag{jj,ii}./shock_min_F);
%
%
%             plot(timeplot,oom_above_min,'k.')
%             hold on
%             title('oom above min force')
%
             %keyboard
            end

        end


        %keyboard
    end
    %title(num2str(evalval))



end
%keyboard





%% force magnitude

% sgtitle(['Max force and velocity over time for each depth ' channelstring 'm channel - ' grainstring ' cm fill' ])
%  legend('show','Location','Best')



% build derivative of force vs channel position using finite difference approximation
Fmag_fdif = zeros(znumgrids,ynumgrids);
pderiv_smooth = zeros(znumgrids,ynumgrids);
pderiv_smooth_deriv = zeros(znumgrids,ynumgrids);



for jj = 1:znumgrids
    % assign depth values to temporary variable u
%     u = Fmag_shocktrack(jj,:);
%     dx = wspacing;
%
%     % use a forward difference method on LHS (first order)
%     Fmag_fdif(jj,1) = (u(2)-u(1))/dx;
%
%     % use backward difference on the RHS (first order)
%     Fmag_fdif(jj,end) = (u(end) - u(end-1))/dx;
%
%     % use central difference (second order) on RHS - 1 and LHS + 1
%     Fmag_fdif(jj,2) = (u(3) - u(1))/(2*dx);
%     Fmag_fdif(jj,end-1) = (u(end) - u(end-2))/(2*dx);
%
%     % use central difference (fourth order) on interior points
%     Fmag_fdif(jj,3:(end-2)) = (-1*u(5:end) + 8*u(4:(end-1)) - 8*u(2:(end-3)) + u(1:(end-4)))/(12*dx);
%
    Fmag_fdif(jj,:) = finite_diff_d1_fun(Fmag_shocktrack(jj,:),wspacing);

    % derivative
%     figure
%     plot(channel,abs(Fmag_fdif(jj,:)))
%     grid on
%     hold on
%     xlabel('Channel')
%     ylabel('Magnitude of dF/dy [N]')
%     set(gca,'Yscale','log')
%     title(['Finite diff approximation to max force derivative vs channel position - ' num2str(jj)])
    % percent of max derivative value
    pderiv = abs(Fmag_fdif(jj,:))/max(abs(Fmag_fdif(jj,:)));
    % smoothing
    pderiv_a = smoothdata(smoothdata(smoothdata(smoothdata(pderiv,'movmean'),'movmean'),'movmean'),'movmean');

    %pderiv_b = smoothdata(smoothdata(pderiv,'movmedian'),'movmedian');
    %pderiv_c = smoothdata(smoothdata(smoothdata(pderiv,'gaussian'),'gaussian'),'movmean');
    %pderiv_d = smoothdata(smoothdata(pderiv,'lowess'),'lowess');
    %pderiv_e = smoothdata(smoothdata(pderiv,'loess'),'loess');
    %pderiv_f = smoothdata(smoothdata(pderiv,'rlowess'),'rlowess');
    %pderiv_g = smoothdata(smoothdata(pderiv,'rloess'),'rloess');
    %pderiv_h = smoothdata(smoothdata(pderiv,'sgolay'),'sgolay');

    pderiv_smooth(jj,:) = pderiv_a; % store
    pderiv_smooth_deriv(jj,:) = finite_diff_d1_fun(pderiv_a,wspacing);
%     % curve fit
%     [curve_p1,gof_p1] = fit(channel.',pderiv.','poly1'); % linear polynomial
%     %[curve_p2,gof_p2] = fit(channel.',pderiv.','poly2'); % quadratic polynomial
%     [curve_l1,gof_l1] = fit(channel.',pderiv.','linearinterp'); % piecewise linear interpolation
%     [curve_c1,gof_c1] = fit(channel.',pderiv.','cubicinterp'); % piecewise cubic interp
%     [curve_s1,gof_s1] = fit(channel.',pderiv.','smoothingspline'); % smoothing spline

end

%


% stand alone plot for looking at the max forces referenced to the steady
% state force
% first array is referenced all to the same steady state force (lowest
% experienced in sim), second is to the steady state force by depth
% keyboard
%Fmag_shocktrack_logref_standardref = log10(Fmag_shocktrack./ref);
Fmag_shocktrack_logref_depthref = zeros(znumgrids,ynumgrids);

%
% figure
% axes('NextPlot','add')
%subplot(1,1,1)
for jj = 1:znumgrids
    % Fmag_rowavg_init,Fmag_rowstd_init
%     refatdepth = Fmag_rowavg_init(jj);
%     Fmag_shocktrack_logref_depthref(jj,:) = log10(Fmag_shocktrack_pos(jj,:)./refatdepth);
%    	thisrow_force = (Fmag_shocktrack_pos(jj,:));
%     plot(channel,thisrow_force,'--o','Color',rowcolorfull(jj,:));
%     xlabel('Radial Position [cm]')
%     ylabel('Force [N]')
%     set(gca, 'YScale', 'log')
%     keyboard



%     subplot(1,2,1)
%     depthstr = ['Depth = ' num2str(depthplot(jj)) ' cm'];
%     %yyaxis left
%     hold on
%     plot(channel,Fmag_shocktrack_logref_standardref(jj,:),'--o','DisplayName',[depthstr],'Color',rowcolorfull(jj,:))
%     grid on
%     ylabel('Max |F|, orders of magnitude relative to global steady state force')
%     xlabel('Radial Position [cm]')
%     %set(gca, 'YScale', 'log')
%     hold off
% %     yyaxis right
% %     hold on
% %     plot(channel,Fmag_shocktrack_time(jj,:),'-','DisplayName',[depthstr],'Color',rowcolorfull(jj,:))
% %     hold off
% %     ylabel('Time [s]')
%     title('Max |F|, relative magnitude compared to SS vs channel position')
%     set(gca,'fontsize',14)

    %subplot(1,2,2)
    depthstr = ['Depth = ' num2str(depthplot(jj)) ' cm'];
    %yyaxis left
%    hold on
    % condition for plotting based on analysis region

%     if (jj >= regionstartdepth) && (jj <= regionstopdepth)
%         plot(channel,Fmag_shocktrack_logref_depthref(jj,:),'--o','DisplayName',[depthstr],'Color',rowcolorfull(jj,:))
%         title('1522')
%     end
%    grid on
%    ylabel('Max |F|, orders of magnitude relative to steady state force by depth')
%    xlabel('Radial Position [cm]')
    %set(gca, 'YScale', 'log')
%    hold off
%     yyaxis right
%     hold on
%     plot(channel,Fmag_shocktrack_time(jj,:),'-','DisplayName',[depthstr],'Color',rowcolorfull(jj,:))
%     hold off
%     ylabel('Time [s]')
    %title('Max |F|, relative magnitude compared to SS (by depth) vs channel position')
%    set(gca,'fontsize',14)

end

% sgtitle(['Max force relative order of magnitude ' channelstring 'm channel - ' grainstring ' cm fill' ])
%title(['Max force relative order of magnitude ' channelstring 'm channel - ' grainstring ' cm fill' ])

%legend('show','Location','Best')


% create the average for the comparison plot
Fmag_logref_mean = nanmean(Fmag_shocktrack_logref_depthref(regionstartdepth:regionstopdepth,:));
Fmag_logref_std = nanstd(Fmag_shocktrack_logref_depthref(regionstartdepth:regionstopdepth,:));


F_shocktrack_mean = nanmean(Fmag_shocktrack_pos(regionstartdepth:regionstopdepth,:));



F_overlaptrack_mean = nanmean(Fmag_shocktrack_maxoverlap(regionstartdepth:regionstopdepth,:));
F_overlaptrack_mean_time = nanmean(Fmag_shocktrack_maxoverlap_time(regionstartdepth:regionstopdepth,:));


F_velocitytrack_mean = nanmean(Fmag_shocktrack_maxvelocity(regionstartdepth:regionstopdepth,:));


% F_overlaptrack_mean_vsdepth = zeros(1,znumgrids);
% F_velocitytrack_mean_vsdepth = zeros(1,znumgrids);
%
%
% for jj = 1:znumgrids
%     F_overlaptrack_mean_vsdepth(1,jj) = nanmean(Fmag_shocktrack_maxoverlap(regionstartdepth:regionstopdepth,jj));
%     F_velocitytrack_mean_vsdepth(1,jj) = nanmean(Fmag_shocktrack_maxvelocity(regionstartdepth:regionstopdepth,jj));
% end




F_overlaptrack_mean_vsdepth = zeros(1,znumgrids);
F_overlaptrack_mean_first_vsdepth = zeros(1,znumgrids);
F_velocitytrack_mean_vsdepth = zeros(1,znumgrids);
F_forcetrack_mean_vsdepth = zeros(1,znumgrids);

% figure
%

% keyboard
% pc = turbo(znumgrids);
for jj = 3:znumgrids-2
    F_overlaptrack_mean_vsdepth(1,jj) = nanmean(Fmag_shocktrack_maxoverlap(regionstartdepth:regionstopdepth,jj));
    F_overlaptrack_mean_first_vsdepth(1,jj) = nanmean(Fmag_shocktrack_maxoverlap_first(regionstartdepth:regionstopdepth,jj));
    F_velocitytrack_mean_vsdepth(1,jj) = nanmean(Fmag_shocktrack_maxvelocity(regionstartdepth:regionstopdepth,jj));
    F_forcetrack_mean_vsdepth(1,jj) = nanmean(Fmag_shocktrack_pos(regionstartdepth:regionstopdepth,jj));

%     plot(abs(Fmag_shocktrack_time(jj,:)),Fmag_shocktrack(jj,:),'.','LineWidth',2','DisplayName',num2str(jj),'Color',pc(jj,:))
%     hold on
%     grid on
%     %set(gca,'Ydir','reverse')
%     ylabel('F [N]')
%     xlabel('t [s]')
%     %title(['Depth vs Wave speed - ' comptype])
%     %title('Propagating c_w')
%     legend('show','Location','best')
%     set(gca,'fontsize',14)
%     set(gca,'yscale','log')

end

    F_mag_avg = mean(Fmag_shocktrack(3:znumgrids-2,:));
    F_mag_avg_time = mean(abs(Fmag_shocktrack_time(3:znumgrids-2,:)));

%F_shocktrack_mean = max(Fmag_shocktrack_pos(regionstartdepth:regionstopdepth,:));

%% wave plotting - using full averages

% ******** radial *********
% ~~~~~~~~~~~ each depth ~~~~~~~~~~~ %
rowcolorfull = turbo(znumgrids);
colcolorfull = turbo(ynumgrids);
Lshort = spmtind + 4; % add some padding so plot doesn't end at min

% ** sanity check below for the averages. turned off for now ** %

% if plotsonind
% figure
% axes('NextPlot','add')
% subplot(1,2,1)
% for j = 1:znumgrids
%     %timeplotradialdepths(j,:) = timeplot.';
%
%     depthstr = ['Depth = ' num2str(depthplot(j)) ' cm'];
%
%     subplot(1,2,1)
%     errorbar(timeplot(1:Lshort),Fmag_rowavg_full(j,1:Lshort),Fmag_rowstd_full(j,1:Lshort),'--o','DisplayName',[depthstr],'Color',rowcolorfull(j,:))
%     hold on
%     grid on
%     set(gca, 'YScale', 'log')
%     xlabel('Time [s]')
%     ylabel('Force [log(N)]')
%     title('Force vs time for each depth')
%
%     subplot(1,2,2)
%     errorbar(timeplot(1:Lshort),100*vy_rowavg_full(j,1:Lshort), 100*vy_rowstd_full(j,1:Lshort), '--o','DisplayName',[depthstr],'Color',rowcolorfull(j,:))
%     hold on
%     grid on
%     xlabel('Time [s]')
%     ylabel('Average Particle Velocity [cm/s]')
%     title('Average Particle Velocity vs time')
%
%
%
% end
%
% sgtitle(['Shock period trends in force and velocity over time for each depth'  channelstring 'm channel - ' grainstring ' cm fill' ])
% legend('show','Location','Best')
%
%
%
% figure
% axes('NextPlot','add')
% subplot(1,2,1)
% for j = 1:ynumgrids
%     %timeplotradialdepths(j,:) = timeplot.';
%
%     channelstr = ['Y = ' num2str(channel(j)) ' cm'];
%
%     subplot(1,2,1)
%     errorbar(timeplot(1:Lshort),Fmag_colavg_full(j,1:Lshort),Fmag_colstd_full(j,1:Lshort),'--o','DisplayName',[channelstr],'Color',colcolorfull(j,:))
%     hold on
%     grid on
%     set(gca, 'YScale', 'log')
%     xlabel('Time [s]')
%     ylabel('Force [log(N)]')
%     title('Force vs time for each radial')
%
%     subplot(1,2,2)
%     errorbar(timeplot(1:Lshort),100*vy_colavg_full(j,1:Lshort), 100*vy_colstd_full(j,1:Lshort), '--o','DisplayName',[channelstr],'Color',colcolorfull(j,:))
%     hold on
%     grid on
%     xlabel('Time [s]')
%     ylabel('Average Particle Velocity [cm/s]')
%     title('Average Particle Velocity vs time')
%
%
%
% end
%
% sgtitle(['Shock period trends in force and velocity over time for each radial -'  channelstring 'm channel - ' grainstring ' cm fill' ])
% legend('show','Location','Best')
%
%
% end
% ~~~~~~~~~~~ maximums ~~~~~~~~~~~ %
% ** force **


if plotsonind
% channel length (force radial)
% h_max_averages = figure;
% subplot(1,2,1)
% % semilogy(max_forcedist_radial_dw,max_force_radial_dw,'--o','LineWidth',2)
% errorbar(max_Fmagdist_radial_fullw,max_Fmag_radial_fullw,max_Fmagstd_radial_fullw,':o','LineWidth',2)
% grid on
% hold on
% errorbar(max_Fxdist_radial_fullw,max_Fxradial_fullw,max_Fxstd_radial_fullw,':o','LineWidth',2)
% errorbar(max_Fydist_radial_fullw,max_Fyradial_fullw,max_Fystd_radial_fullw,':o','LineWidth',2)
% errorbar(max_Fzdist_radial_fullw,max_Fzradial_fullw,max_Fzstd_radial_fullw,':o','LineWidth',2)
% set(gca, 'YScale', 'log')
% xlabel('Y [m]')
% ylabel('Force [N]')
% hold off
% yyaxis right
% plot(max_Fmagdist_radial_fullw,max_Fmagtime_radial_fullw)
% ylabel('Time [s]')
% title('Depth averaged maximum force and time vs channel position')
% legend('Force (magnitude)','Force (X)','Force (Y)','Force (Z)', 'Time of max force','Location','Best')
% set(gca,'fontsize',14)
%
%
% % ** velocity **
%
%
% % channel length (velocity radial)
% subplot(1,2,2)
% errorbar(max_vmagdist_radial_fullw,100*max_vmag_radial_fullw,100*max_vmagstd_radial_fullw,'--o','LineWidth',2)
% hold on
% grid on
% errorbar(max_vxdist_radial_fullw,100*max_vx_radial_fullw,100*max_vxstd_radial_fullw,'--o','LineWidth',2)
% errorbar(max_vydist_radial_fullw,100*max_vy_radial_fullw,100*max_vystd_radial_fullw,'--o','LineWidth',2)
% errorbar(max_vzdist_radial_fullw,100*max_vz_radial_fullw,100*max_vzstd_radial_fullw,'--o','LineWidth',2)
% %plot(max_vydist_radial_fullw,100*max_vy_radial_fullw,':o','LineWidth',2)
% xlabel('Y [m]')
% ylabel('Velocity [cm/s]')
% hold off
% yyaxis right
% plot(max_vydist_radial_fullw, max_vytime_radial_fullw)
% title('Depth averaged maximum velocity vs channel position')
% legend('Velocity (magnitude)','Velocity (X)', 'Velocity (Y)', 'Velocity (Z)','Time of max velocity', 'All windows','Location','Best')
% set(gca,'fontsize',14)
%
% set(h_max_averages,'Units', 'pixels', 'Position', [500, 200, 1200, 600]);
%
% sgtitle(['Maximums over time - '  channelstring 'm channel - ' grainstring ' cm fill'])

end



%% find the runout distance from the intersection method (new 04/26/2021)

% find the index that the shock mintime occurs at
runoutindex_F = zeros(znumgrids,1);
shock_runout_F = zeros(znumgrids,1);
shock_speed_F = zeros(znumgrids,1);
shockmin_F = zeros(znumgrids,1);
shockminstd_F = zeros(znumgrids,1);

runoutindex_V = zeros(znumgrids,1);
shock_runout_V = zeros(znumgrids,1);
shock_speed_V = zeros(znumgrids,1);


shockspeedslope = zeros(znumgrids,1);% 1/slopeint(2); % m/s
shockspeedslope_exp = zeros(znumgrids,1);% 1/slopeint(2); % m/s
shocktimecalc = zeros(znumgrids,1);
shockspeedslope_arrival = zeros(znumgrids,1);
shockspeedslope_maxforce = zeros(znumgrids,1);

lastforceinds = zeros(znumgrids,1);
fitparamstore = zeros(znumgrids,3);

wave_measurement = 0;
wave_time = 0;
wave_fit = 0;
wave_fit_time = 0;


for jj = 1:znumgrids
    if (jj > regionstartdepth) && (jj <= regionstopdepth)
%         lastforceindstemp = find(abs(pderiv_smooth_deriv(jj,:)) <= (10^-4),1,'first');
        % call function to that finds runout distance
        [terminationp,lastforceindstemp,p] = terminationfun(channel,Fmag_shocktrack_pos(jj,:),wspacing,jj,Fmag_rowavg_init(jj));
        % store fit parameters
        fitparamstore(jj,:) = p;



        if isempty(lastforceindstemp) || lastforceindstemp == 0

            lastforceinds(jj,1) = NaN;
            runoutindex_F(jj,1) = NaN;
            shock_runout_F(jj,1) = NaN; % m
            shocktimecalc(jj,1) = NaN;
            shock_speed_F(jj,1) = NaN;
            shockmin_F(jj,1) = NaN;
            shockminstd_F(jj,1) = NaN;
            shockspeedslope(jj,1) = NaN;
            shockspeedslope_exp(jj,1) = NaN;
            shockspeedslope_arrival(jj,1) = NaN;

        else
            lastforceinds(jj,1) = lastforceindstemp;
            runoutindex_F(jj,1) = lastforceindstemp;
            shock_runout_F(jj,1) = channel(lastforceindstemp)/100; % m
            % wave speed 1 calc - use only the first 50 cm
            %lastforceindstemp = shockstop;

            if jj == 5
                %keyboard
            end

            if shocktrigger(jj,1) < shockstop
                stopcalc = shocktrigger(jj,1);
            else
                stopcalc = shockstop;
            end
            %stopcalc = shocktrigger(jj,1)
            stopcalc = shockstop;

            % exclude any data with a max at time = 0

            ydata = Fmag_shocktrack_time(jj,1:stopcalc);
            xdata = channel(1:stopcalc)./100; % convert back to m

            %keyboard

            range_prop = 40; %80;
            range_exp = 35;
            range_SW = 20; %35;
            % reduce region for shocks
            if vel > 1
                range_exp = 12; %30;
            end

            if r_p == 5E-3 || r_p == 0.5E-3 || r_p == 125E-6
                range_exp = 12;
            end

            if r_p == 125E-6
                range_exp = 12;
            end

            if r_p == 125E-6
                range_prop = 90;
                range_SW = 60;
            end

            % exclude the row if it didn't get populated (pressurized case)
            if sum(isnan(Fmag_shocktrack_time(jj,2:range_exp))) > 0
                % set this as the exponential decay regime speed

                [time_uniform,fit_uniform,wavespeed,fitparam] = ...     % 20:stopcalc
                    least_fit_fun_nearestneighbor(Fmag_shocktrack_time(jj,2:range_exp),channel(2:range_exp)./100,vel,jj);



                % set this as the down channel speed, not the exp (121282022)
                [time_uniform_exp,fit_uniform_exp,wavespeed_exp,~] = ...
                    least_fit_fun_nearestneighbor(Fmag_shocktrack_time(jj,range_SW:range_prop),channel(range_SW:range_prop)./100,vel,jj);
            else
                time_uniform = nan;
                fit_uniform = nan;
                wavespeed = nan;
                wavespeed_exp = nan;
            end
            %keyboard

            if jj == sz
                wave_measurement = channel(2:stopcalc)./100;
                wave_time = Fmag_shocktrack_time(jj,2:stopcalc);
                wave_time(wave_measurement == 0) = nan;
                wave_measurement(wave_measurement == 0) = nan;
                wave_fit = fit_uniform;
                wave_fit_time = time_uniform;

            end

            %timeinds = find(ydata2 < 0.2);


            X = [ones(length(xdata),1) xdata.'];
            slopeint = X\ydata.'; % time (s)/ m
%             shockspeedslope(jj,1) = 1/slopeint(2);

            shockspeedslope(jj,1) = wavespeed(1);
            shockspeedslope_exp(jj,1) = wavespeed_exp(1);

            %shockspeedslope_expreg(jj,1) =;

%             shocktimecalc(jj,1) = slopeint(1) + slopeint(2)*shock_runout_F(jj,1);
            shocktimecalc(jj,1) = 0;
            %shocktimecalc(jj,1) = Fmag_shocktrack_time(jj,lastforceindstemp);
            shock_speed_F(jj,1) = shockspeedslope(jj,1); %shock_runout_F(jj,1)/shocktimecalc(jj,1);
            % find the index in time array
            [~,shockmindex] = min(abs(timeplot - shocktimecalc(jj,1)));
%             shockmin_F(jj,1) = Fmag_rowavg_full(jj,shockmindex);
%             shockminstd_F(jj,1) = Fmag_rowstd_full(jj,shockmindex);
            shockmin_F(jj,1) = 0;
            shockminstd_F(jj,1) = 0;


            % repeat wave speed for one crossing
            startcalc = 1;
            ydata3 = Fmag_shocktrack_time_onecross(jj,startcalc:stopcalc);
            xdata3 = xdata(startcalc:end);
            X3 = [ones(length(xdata3),1) xdata3.'];
            % remove entries that have the first time index (time = 0)
            non_z_inds = find(ydata3);

            % if non zero entries proceed as normal
            if length(non_z_inds) == length(ydata3)
                slopeint2 = X3\ydata3.'; % time (s)/ m
                shockspeedslope_arrival(jj,1) = 1/slopeint2(2);

%                 figure
%                 plot(ydata3,X(:,2),'k.')
%                 hold on
%                 ylabel('position')
%                 xlabel('time')
%
%                 title(['arrival: v ' num2str(evalval) ' depth index ' num2str(jj)])
%                 keyboard
            else % remove the zero entry
                ydata3_new = ydata3(non_z_inds);
                xdata3_new = xdata3(non_z_inds);
                X3_new = [ones(length(xdata3_new),1) xdata3_new.'];

                slopeint2 = X3_new\ydata3_new.'; % time (s)/ m
                shockspeedslope_arrival(jj,1) = 1/slopeint2(2);

%                 figure
%                 plot(ydata3_new,X3_new(:,2),'k.')
%                 hold on
%                 ylabel('position')
%                 xlabel('time')
%
%                 title(['arrival: v ' num2str(evalval) ' depth index ' num2str(jj)])
%                 keyboard

            end



           % repeat for max force
%             figure
            ydata4 = Fmag_shocktrack_time_maxforce(jj,startcalc:stopcalc);
            slopeint3 = X3\ydata4.'; % time (s)/ m
            shockspeedslope_maxforce(jj,1) = 1/slopeint3(2);
%
%             keyboard
            if jj == 12
                %keyboard
            end

%             figure
%             plot(ydata4,X3(:,2),'k.')
%             hold on
%             ylabel('position')
%             xlabel('time')
%
%             title(['max: v ' num2str(evalval) ' depth index ' num2str(jj)])
%         keyboard
%         % for paper - plotting of time calculation
%         [~,~,~] = terminationfun2(channel,Fmag_shocktrack(jj,:),wspacing,jj);
%         h_termination_time = figure;
        %set(0, 'CurrentFigure', h_termination)
        if jj == 5
%             h_channelvstime = figure;
%             subplot(1,2,2)
%             plot(Fmag_shocktrack_time(jj,:),channel./100,'k.','LineWidth',2,'DisplayName',['Time of Max Force at depth ' num2str(jj) ' cm'])
%             hold on
%             grid on
%             plot((slopeint(1) + slopeint(2)*xdata),xdata,'b--','LineWidth',2,'DisplayName','Fit')
%             plot(shocktimecalc(jj,1),channel(lastforceindstemp)/100,'r*','LineWidth',3,'MarkerSize',10,'DisplayName','Termination')
%             ylabel('Channel Position [m]')
%             xlabel('Time [s]')
%             %title(['Time of max force occurence vs channel position for depth ' num2str(round(depthplot(jj))) 'cm'])
%             legend('show','location','best')
%             set(gca,'fontsize',14)
%             keyboard
        end
%
%         keyboard


%         %plotpos = get(h_termination,'position');
%         % double the width of the plot
%         %set(h_termination,'Units', 'pixels', 'Position', [plotpos(1:2), plotpos(3)*2, plotpos(4)]);
%         keyboard
        %saveas(h_termination,[plotpath 'terminationtimeanddfdy.png']);
%         saveas(h_termination_pos,[plotpath 'termination_pos.png'])
%         saveas(h_termination_time,[plotpath 'termination_time.png'])
        end % end the is empty check
    end % end conditional for in analysis range
end % end over all depths
% set any zero entries in arrival to NaN
shockspeedslope_arrival(shockspeedslope_arrival == 0) = NaN;
shockspeedslope_maxforce(shockspeedslope_maxforce == 0) = NaN;

%keyboard
% set zeros to nan for plotting
shock_runout_F(shock_runout_F==0)=NaN;
shockspeedslope(shockspeedslope==0)=NaN;
shockspeedslope_exp(shockspeedslope_exp==0)=NaN;
shockmin_F(shockmin_F==0)=NaN;
shockminstd_F(shockminstd_F==0)=NaN;
fitparamstore(fitparamstore==0)=NaN;

parammean = nanmean(fitparamstore);
paramstd = nanstd(fitparamstore);



zreg = zeros(1,length(regionstartdepth:regionstopdepth));
%% plot the shockstuffs times and the max mins
if 0


h_shockdetailsF = figure;
subplot(1,3,1)
errorbar(shockmax(regionstartdepth:regionstopdepth),depthplot(regionstartdepth:regionstopdepth),zreg,zreg,shockmaxstd(regionstartdepth:regionstopdepth),shockmaxstd(regionstartdepth:regionstopdepth),'--o','LineWidth',2)
hold on
grid on
errorbar(shockmin_F(regionstartdepth:regionstopdepth),depthplot(regionstartdepth:regionstopdepth),zreg,zreg,shockminstd_F(regionstartdepth:regionstopdepth),shockminstd_F(regionstartdepth:regionstopdepth),'--x','LineWidth',2)
set(gca, 'XScale', 'log')
set(gca,'Ydir','reverse')
xlabel('Force [log(N)]')
ylabel('Depth [cm]')
title('Depth vs Shock max and min')
set(gca,'fontsize',14)


subplot(1,3,2)
plot(shock_runout_F,depthplot,'--o','LineWidth',2)
hold on
grid on
ylabel('Depth [cm]')
xlabel('Radial distance [m]')
set(gca,'Ydir','reverse')
title('Depth vs Shock run out distance')
set(gca,'fontsize',14)

subplot(1,3,3)
%plot(depthplot,shock_speed_F,'--o','LineWidth',2)
plot(shockspeedslope,depthplot,'--o','LineWidth',2)
hold on
grid on
% plot(depthplot,shockspeedslope,'--x','LineWidth',2)
ylabel('Depth [cm]')
xlabel('Velocity [m/s]')
set(gca,'Ydir','reverse')
title('Depth vs Shock speed')
set(gca,'fontsize',14)

sgtitle(['Depth vs, (F) Shock strength, duration, speed and runout- ' channelstring 'm channel - ' grainstring ' cm fill' ])

set(h_shockdetailsF,'Units', 'pixels', 'Position', [500, 200, 1250, 450]);

end


%% bulk dilation

% find crater properties
[fc_row,fc_col] = find(ssdilationind_full == 150);
% find first column where crater no longer exists (last consecutive
% void
% fc_col_u = unique(fc_col);
% fc_col_u_binary = [fc_col_u(1); fc_col_u(2:end) - fc_col_u(1:end-1)];
% % find the first non consecutive index
% r_crater_ind = fc_col_u(find(fc_col_u_binary > 1,1,'first') - 1);
%     if isempty(r_crater_ind)
%         r_crater_ind = fc_col_u(end) + 1;
%     end
% r_crater = channel(r_crater_ind) %cm total crater radius
% r_crater_inrange = channel(r_crater_ind)-channel(1); %cm
% % find the lowest row for an angle measurement
% fc_row_one = fc_row(find(fc_col == fc_col(1)));
% z_crater_ind = max(fc_row_one) + 1;
% % z_crater = depthplot(z_crater_ind); % cm
% % % angle
% %
% % craterslope = atand(z_crater/r_crater_inrange);

r_crater_ind = 5;

r_bulk_start = 21; % 10 cm
r_bulk_stop = znumgrids - 50;%241; % 1.2 m (120 cm)

% calculate bulk dilation average for each depth with [rcrater,shockrunout]
bulk_dilation_vsdepth = zeros(znumgrids,1);
bulk_dilation_stddv_vsdepth = zeros(znumgrids,1);

bulk_dilation_vsradial = zeros(ynumgrids,1);
bulk_dilation_stddv_vsradial = zeros(ynumgrids,1);

for jj = 1:znumgrids

    if isnan(lastforceinds(jj))
        bulk_dilation_vsdepth(jj,1) = NaN;
        bulk_dilation_stddv_vsdepth(jj,1) = NaN;
    else
%         bulk_dilation_vsdepth(jj,1) = nanmean(dilationpdiff_SS(jj,r_crater_ind:floor(lastforceinds(jj)/8)));
%         bulk_dilation_stddv_vsdepth(jj,1) = nanstd(dilationpdiff_SS(jj,r_crater_ind:floor(lastforceinds(jj)/8)));
%         bulk_dilation_vsdepth(jj,1) = nanmean(dilationpdiff_SS(jj,r_bulk_start:r_bulk_stop));
%         bulk_dilation_stddv_vsdepth(jj,1) = nanstd(dilationpdiff_SS(jj,r_bulk_start:r_bulk_stop));
    end

end

for jj = 1:ynumgrids
%      %bulk_dilation_vsradial(jj,1) = nanmean(dilationpdiff_SS(:,jj));
%      bulk_dilation_vsradial(jj,1) = sum(dilationpdiff_SS(:,jj));
%      bulk_dilation_stddv_vsradial(jj,1) = nanstd(dilationpdiff_SS(:,jj));
end

%% "core" samples
% plot dilation vs depth from a "core sample" at some selection of depths

% add the expected trend to the dilation plot
rhod = 1900; % kgm-3, density at infinite depth
rho0 = 900; % kgm-3, density at surface
H1 = 5.7; % cm
H2 = 8.6; % cm

z = 0:.1:40; % top 30 cm
L = length(z);
deltarho = zeros(1,L);
d_background = zeros(1,L);
d_coldspot = zeros(1,L);

for jj = 1:L
    [deltarho(1,jj), d_background(1,jj), d_coldspot(1,jj)] = ...
        densityVdepth(rhod,rho0,H1,H2,z(jj));
end
%
% h_depth_core = figure;
% axes('NextPlot','add');
% colrcount = 1;
% corerange = 1:5:100; %1:5:50
% plotcolr = turbo(length(corerange));
%
% for jj = corerange
%
%     templotdil = dilationpdiff_SS(:,jj);
%
%     plot(templotdil(templotdil<0),depthplot(templotdil<0),'o','LineWidth',2,'DisplayName',['R = ' num2str(channel(jj)) ' m'],'Color',plotcolr(colrcount,:))
%     hold on
%     grid on
%     %errorbar(ss_compact_row_avg_yval,ss_compact_row_avg,ss_compact_row_std,'--o','LineWidth',2)
%     ylabel('Depth [cm]')
%     xlabel('Percent Dilation')
%     set(gca,'Ydir','reverse')
%     title(['Depth vs Dilation - ' comptype ' - ' evalval])
%     set(gca,'fontsize',14)
%
%     colrcount = colrcount + 1;
% end
%
% %set(0, 'CurrentFigure', h_dil_vsdepth)
% plot(100*deltarho,z,'k-','LineWidth',2,'DisplayName','Expected')
% legend('show','location','best')
% xlim([-10, 2])
%
% saveas(h_depth_core, [plotpath 'dilationvsdepth_core.png'])

%% Force vs time plotting for specified bins
% set bin location
bd = 15; %15cm in depth
br = 30; %30cm from impact wall
% call fitting function and pass the bin's F vs t

% [~] = dispersionfun(store_Fmag{bd,20},timeplot,channel(20),depthplot(bd));
% [~] = dispersionfun(store_Fmag{bd,br},timeplot,channel(br),depthplot(bd));
% [~] = dispersionfun(store_Fmag{bd,40},timeplot,channel(40),depthplot(bd));


% axes('NextPlot','add')
% subplot(1,1,1)
% force magnitude plot
% for jj = 1:znumgrids
%     figure
%     plot(timeplot,store_Fmag{jj,400},'--o','LineWidth',2)
%     hold on
%     grid on
%     % plot(depthplot,shockspeedslope,'--x','LineWidth',2)
%     ylabel('Force [N]')
%     xlabel('Time [s]')
%     set(gca,'yscale','log')
%     set(gca,'fontsize',14)
%     legend('show','location','best')
%     title(['Depth = ' num2str(jj)])
% end
%
%  keyboard
%% Heat Plots

% use  already created radial distance and depth arrays for x and y labels

% depthstart = ;
% depthstop = ;
lengthstart = 1;
depthstarthm = 3;
% lengthstop = ;

% dilation indicator

%keyboard

% h_dilind = figure;
% % heatmap(100*max_vydist_radial_fullw,100*max_vydepth_fullw,ssdilationind_full)
% imagesc(max_vydist_radial_fullw,100*max_vydepth_fullw,ssdilationind_full)
% xlabel('Radial Distance [m]')
% ylabel('Depth [cm]')
% %title(['Steady State Dilation Indicator - '  channelstring 'm channel - ' grainstring ' cm fill' ])
% colormap(bone)
% colorbar('Ticks',[-100, -20, 50, 150],...
%          'TickLabels',{'Compaction','Partial Dilation','Dilation','Void'})
% set(gca,'fontsize',14)
% set(h_dilind,'Units', 'pixels', 'Position', [500, 200, 1000, 400]);

% save figure
%saveas(h_dilind, [plotpath 'dilationindicator.png'])



%% heights
% single plot
% what does a single height vs time look like at some R?
R_height_check = 100;
% h_height_vs_time = figure;
% axes('NextPlot','add');


% for jj = 10:ynumgrids
%    for kk = 1:Ltp
%
%
%    axes('NextPlot','add');
%    subplot(1,1,1)
%    plot(timeplot(kk),store_heights{1,jj}(kk),'k.');
%    %plot(timeplot(kk),store_heights{1,R_height_check}(kk),'k.');
%    hold on
%    grid on
%    xlabel('Time [s]')
%    ylabel('z [cm]')
%    set(gca,'FontSize',14)
%    end
%    %keyboard
%
% end


store_max_height = zeros(1,ynumgrids);

for jj = 1:Ltp
    tstep_heights = store_heights{1,jj};
    check = tstep_heights > store_max_height;
    store_max_height(check) = store_heights{1,jj}(check);
end
% take the pdiff later after the init height is stored


% multi plot
% what does a single height vs time look like at many R?
% also, plot the initial and final heights
store_init_height = zeros(1,ynumgrids);
store_final_height_mean = zeros(1,ynumgrids);
store_final_height_std = zeros(1,ynumgrids);

store_heights_bin = zeros(1,Ltp);

for ii = 1:ynumgrids

%print height vs time for every N radial
% plot for every N radials
if ii == sr %   ~rem(ii,50)
%     figure
%     axes('NextPlot','add');
    for kk = 1:Ltp
       %axes('NextPlot','add');
%        subplot(1,1,1)
%        plot(timeplot(kk),store_heights{1,kk}(ii),'k.');
%        hold on
%        grid on
%        xlabel('Time [s]')
%        ylabel('z [m]')
%        set(gca,'FontSize',14)
       %title(['R = ' num2str(channel(ii)) ' cm'])

       % for output plotting of height vs time comparison plot
       store_heights_bin(1,kk) = store_heights{1,kk}(ii);

    end
    % keyboard
end


% set(0, 'CurrentFigure', h_height_BegandEnd)
% axes('NextPlot','add')
% subplot(1,1,1)
% plot(channel(ii),

% just store the initial and final heights here, then plot after

store_init_height(1,ii) = store_heights{1,1}(ii);
temp_height = zeros(1,length((Ltp-5):Ltp)); % 20

for rr = (Ltp-5):Ltp % 20
    temp_height(1,rr - (Ltp-6)) = store_heights{1,rr}(ii); % 21
end

store_final_height_mean(1,ii) = mean(temp_height);
store_final_height_std(1,ii) = std(temp_height);

%keyboard
%store_max_height(1,ii) = max(store_heights{1,ii});
%
% if ii == 108
%     keyboard
% end


end

% pdiff of max height
store_max_height_pdiff = 100*(store_max_height - store_init_height)./store_init_height;

% h_height_BegandEnd = figure;
% plot(channel/100,store_init_height,'ko','DisplayName','Initial Height')
% hold on
% grid on
% plot(channel/100,store_final_height_mean,'bo','DisplayName','Final Height')
% ylabel('z [m]')
% xlabel('R [m]')
% legend('show','location','best')
% title(['Bed Height - ' evalval])

%keyboard
height_pdiff_out = 100*(store_final_height_mean-store_init_height)./store_init_height;
height_out = store_final_height_mean-store_init_height;
% keyboard
% h_height_pdiff = figure;
% plot(channel/100,height_pdiff_out,'ko','DisplayName','Initial Height')
% hold on
% grid on
% ylabel('% difference z')
% xlabel('R [m]')
% legend('show','location','best')
% title(['Bed Height % difference - ' evalval])



% h_heightvariance_BegandEnd_ = figure;
% plot(channel/100,10^3 *store_final_height_std,'ko','DisplayName','\sigma')
% hold on
% grid on
% ylabel('10^3 \sigma z [m]')
% xlabel('R [m]')
% legend('show','location','best')
% title(['Bed Height std - ' evalval])

% loop for steady state height and variance for each cell

grid_ss_height_mean = zeros(znumgrids,ynumgrids);
grid_ss_height_std = zeros(znumgrids,ynumgrids);
grid_init_height = zeros(znumgrids,ynumgrids);
grid_ss_height_pdiff = zeros(znumgrids,ynumgrids);
%h_depth_percent = figure

depthcolor = turbo(znumgrids);
% h_heightpercentchange_depth = figure;
%
% for jj = 1:znumgrids
%     for kk = 1:ynumgrids
%
%         grid_init_height(jj,kk) = store_psensor_z{jj,kk}(1);
%         temp_heights = [];
%         temp_heights = store_psensor_z{jj,kk}((Ltp-50):Ltp);
%
%         grid_ss_height_mean(jj,kk) = mean(temp_heights);
%         grid_ss_height_std(jj,kk) = std(temp_heights);
%
%         grid_ss_height_pdiff(jj,kk) =  100*(grid_ss_height_mean(jj,kk) - grid_init_height(jj,kk)) ...
%             / grid_init_height(jj,kk);
%
%     end % end y grids
%
%     plot(channel/100,grid_ss_height_pdiff(jj,:),'k.--','Color',depthcolor(jj,:));
%     grid on
%     hold on
%
% end % end z grids
% xlabel('R [m]')
% ylabel('\Delta Z [%]')
% set(gca,'fontsize',14')
% title(num2str(evalval))
%keyboard

% store_final_height_mean(1,ii) = mean(temp_height);
% store_final_height_std(1,ii) = std(temp_height);

%keyboard
% steady state dilation
% figure
% heatmap(100*max_vydist_radial_fullw(lengthstart:lengthstop),100*max_vydepth_fullw(depthstarthm:depthstop),dilationpdiff_SS(depthstarthm:depthstop,lengthstart:lengthstop))
% xlabel('Radial Distance [cm]')
% ylabel('Depth [cm]')
% title('Steady State Dilations')

% max dilation
% figure
% heatmap(100*max_vydist_radial_fullw(lengthstart:lengthstop),100*max_vydepth_fullw(depthstarthm:depthstop),dilationpdiff_max(depthstarthm:depthstop,lengthstart:lengthstop))
% xlabel('Radial Distance [cm]')
% ylabel('Depth [cm]')
% title('Maximum Dilation Experienced')



datetime.setDefaultFormats('default','yyyy-MM-dd')
currenttime = datetime;
formatOut = 'mmddyyyy';
timestring = datestr(currenttime,formatOut);
grainstring = num2str(round(100*grainline));
channelstring = num2str(channellength);
porstring = num2str(round(100*porosity));

%% Particle sensor plots

if 0 % on off control
    ps_color = turbo(znumgrids);
    h_psensor_force = figure;
    axes('NextPlot','add');
    h_psensor_normalstressmag = figure;
    axes('NextPlot','add');
    h_psensor_shearstressmag = figure;
    axes('NextPlot','add');
    h_psensor_vel_mag = figure;
    axes('NextPlot','add');
    h_psensor_stressyy = figure;
    axes('NextPlot','add');
    h_psensor_stressyz = figure;
    axes('NextPlot','add');
    h_psensor_positionz = figure;
    axes('NextPlot','add');



    for kk = 1:znumgrids

       psensor_fmag = (store_psensor_Fx{kk,1}.*store_psensor_Fx{kk,1} + ...
           store_psensor_Fz{kk,1}.*store_psensor_Fz{kk,1} + ...
           store_psensor_Fz{kk,1}.*store_psensor_Fz{kk,1}).^0.5;

       set(0, 'CurrentFigure', h_psensor_force)
       axes('NextPlot','add')
       subplot(1,1,1)
       semilogy(timeplot,psensor_fmag,'o','DisplayName',num2str(store_psensor_z{kk,1}(1)),'Color',ps_color(kk,:));
       hold on
       grid on
       xlabel('Time [s]')
       ylabel('Force [N]')
       set(gca,'FontSize',14)
       set(gca, 'YScale', 'log')



       psensor_normal_stress_mag = (store_psensor_stress_xx{kk,1}.*store_psensor_stress_xx{kk,1} + ...
           store_psensor_stress_yy{kk,1}.*store_psensor_stress_yy{kk,1} + ...
           store_psensor_stress_zz{kk,1}.*store_psensor_stress_zz{kk,1}).^0.5;

       set(0, 'CurrentFigure', h_psensor_normalstressmag)
       axes('NextPlot','add')
       subplot(1,1,1)
       semilogy(timeplot,psensor_normal_stress_mag-psensor_normal_stress_mag(1),'o','DisplayName',num2str(store_psensor_z{kk,1}(1)),'Color',ps_color(kk,:));
       hold on
       grid on
       xlabel('Time [s]')
       ylabel('Normal stress [N*m]')
       set(gca,'FontSize',14)
       set(gca, 'YScale', 'log')


       psensor_shear_stress_mag = (store_psensor_stress_xy{kk,1}.*store_psensor_stress_xy{kk,1} + ...
           store_psensor_stress_xz{kk,1}.*store_psensor_stress_xz{kk,1} + ...
           store_psensor_stress_yz{kk,1}.*store_psensor_stress_yz{kk,1}).^0.5;

       set(0, 'CurrentFigure', h_psensor_shearstressmag)
       axes('NextPlot','add')
       subplot(1,1,1)
       semilogy(timeplot,psensor_shear_stress_mag-psensor_shear_stress_mag(1),'o','DisplayName',num2str(store_psensor_z{kk,1}(1)),'Color',ps_color(kk,:));
       hold on
       grid on
       xlabel('Time [s]')
       ylabel('Shear stress [N*m]')
       set(gca,'FontSize',14)
       set(gca, 'YScale', 'log')


      psensor_vel_mag = (store_psensor_VX{kk,1}.*store_psensor_VX{kk,1} + ...
           store_psensor_VY{kk,1}.*store_psensor_VY{kk,1} + ...
           store_psensor_VZ{kk,1}.*store_psensor_VZ{kk,1}).^0.5;

       set(0, 'CurrentFigure', h_psensor_vel_mag)
       axes('NextPlot','add')
       subplot(1,1,1)
       plot(timeplot,psensor_vel_mag,'o','DisplayName',num2str(store_psensor_z{kk,1}(1)),'Color',ps_color(kk,:));
       hold on
       grid on
       xlabel('Time [s]')
       ylabel('Velocity [m/s]')
       set(gca,'FontSize',14)



       set(0, 'CurrentFigure', h_psensor_stressyy)
       axes('NextPlot','add')
       subplot(1,1,1)
       semilogy(timeplot,abs(store_psensor_stress_yy{kk,1}),'o','DisplayName',num2str(store_psensor_z{kk,1}(1)),'Color',ps_color(kk,:));
       hold on
       grid on
       xlabel('Time [s]')
       ylabel('Normal stress [N*m]')
       set(gca,'FontSize',14)
       set(gca, 'YScale', 'log')


       set(0, 'CurrentFigure', h_psensor_stressyz)
       axes('NextPlot','add')
       subplot(1,1,1)
       plot(timeplot,store_psensor_stress_yz{kk,1},'o','DisplayName',num2str(store_psensor_z{kk,1}(1)),'Color',ps_color(kk,:));
       hold on
       grid on
       xlabel('Time [s]')
       ylabel('Shear stress [N*m]')
       set(gca,'FontSize',14)


       set(0, 'CurrentFigure', h_psensor_positionz)
       axes('NextPlot','add')
       subplot(1,1,1)
       plot(timeplot,100*(store_psensor_z{kk,1}-store_psensor_z{kk,1}(1)),'o','DisplayName',num2str(store_psensor_z{kk,1}(1)),'Color',ps_color(kk,:));
       hold on
       grid on
       xlabel('Time [s]')
       ylabel('Change in height [cm]')
       set(gca,'FontSize',14)

    end

    keyboard
end
%
end

function [drhoz,drho1,drho2] = densityVdepth(rhod,rho0,H1,H2,z)

drho1 = rhod - (rhod - rho0) * exp(-z/H1);
drho2 = rhod - (rhod - rho0) * exp(-z/H2);


drhoz = (drho2)/(drho1) - 1 ;

end
