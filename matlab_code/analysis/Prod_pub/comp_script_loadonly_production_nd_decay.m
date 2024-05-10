% comp_script_load_only_production_nd_decay.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

% This wrapper script manages the analysis of the results of lateral piston impacts
% within a granular channel (simulations performed in LIGGGHTS). This script does not
% operate directly on the LIGGGHTS output files, but instead takes as input a
% pre-processed *.mat file that represents the average particle states over time for
% an array of virtual sensors (bins). For details on the creation of the *.mat file
% or the simulation methodology, see our public repository (link at the top). This
% wrapper script performs further averaging across the virtual bins and produces
% average trends within the channel vs time (particle velocity, particle force,
% particle stresss, bulk packing fraction, wave runout distances). As noted in the
% README, this script contains additional experimental functionality that has been
% been left intact, but commented out (for example, comparison of bed height over time).
% We include this additional material as a reference for those implementing a similar
% methodology. Therefore, I describe only the relevant analysis outputs (those not
% commented out) below.
%
% % ----------- Cases ----------- %
% Three sample cases are included in this work, representing different material parameters
% of the particles in the channel or piston characteristics. The cases are:
% 1. A particle size sweep (changing radius, R) using 'hard' particles (E = 50 GPa,
% rho = 3 kg/m3), 'size_sweep_medium_rho3000_E50GPa_decay'
% 2. A piston mass sweep using 'boulder' particles (R = 1.25 m, E = 50 GPa,
% rho = 3 kg/m3), 'decay_boulder_impactorMass'
% 3. A piston velocity sweep using 'boulder' particles (R = 1.25 m, E = 50 GPa,
% rho = 3 kg/m3), 'decay_boulder_impactorMass'
%
%
% % ----------- Workflow ----------- %
% This code operates on a single case at a time. The infrastructure for the three cases
% above is provided below and can be used as a template for any additional analysis. Only
% one of the header details below should be left uncommented. The steps to producing the
% output are:
% 1. Ensure the header information of the case of interest is uncommented and comment
% out the header detail for all other cases
% 2. Modify any of the header detail as needed (path locations, ordering of files,
% labels, any manually identified runout distances you want to plot)
% 3. Ensure that the case name (comptype variable, below) exists in the ProductionPlots
% direction (ex: ProductionPlots/decay_boulder_velocity)
% 4. Run comp_script_load_only_production_nd_decay
%
%
% % ----------- Output ----------- %
% This wrapper script will produce 9 plots that correspond to some of the figures
% produced in our works as well as a few diagnostic plots. The plots are described
% in detail in the README, 'Code output' section. A 'write_file' which is a text file
% containing the most salient details from the test run is also written to the write_files
% directory, named for the test case considered (and time stamped)
%


%% path
% path to stored *.mat file directory. A link to sample *.mat files used in our work
% is provided in the README
matfile_location = "/Users/hartzell_lab/Documents/COLDSPOTS/MATLAB ANALYSIS/ShockInducedDilation/Prod_pub/mat_files/";


%% size - medium beds - decay - rho 3000 - E 50 GPa

% file location
 matfile_arr =[strcat(matfile_location, "size_sweep_medium_rho3000_E50GPa_decay/2m_channel_20cmfill_medium_1kms_decay_size_1.25mm_rho3000_E50GPa_03072024.mat");
    strcat(matfile_location, "size_sweep_medium_rho3000_E50GPa_decay/20m_channel_2mfill_medium_1kms_decay_size_12.5mm_rho3000_E50GPa_03082024.mat");
     strcat(matfile_location, "size_sweep_medium_rho3000_E50GPa_decay/200m_channel_20mfill_medium_1kms_decay_size125mm_rho3000_E50GPa_02282024.mat");
     strcat(matfile_location, "size_sweep_medium_rho3000_E50GPa_decay/2km_channel_200mfill_medium_1kms_size_1.25m_rho3000_E50GPa_03062024.mat");];

Lf = length(matfile_arr);
% legend entries, corresponding to the investigated parameter
evalvals = ["R = 1.25 mm","R = 12.5mm", "R = 125 mm","R = 1.25 m"];
% piston velocity, m/s
vel = [1000,1000,1000,1000];
% test/case name
comptype = 'size_sweep_medium_rho3000_E50GPa_decay';
% printing flag (vestigial)
velyesno = 0;
% index range over which to measure dilation (vestigial)
r_start_ind = 61;
r_end_buff = 71;
% flag for considering total particle overlap, not just average (vestigial)
tot_deltas_onoff = 0;

% elastic modulii
E_p = [50E9,50E9,50E9,50E9]; % Pa
% radii
r_p = [1.25E-3,12.5E-3,125E-3,1.25 ]; % m
% densities
rho_p = [3000,3000, 3000,3000]; % kg/m3
% particle-particle cohesion (SJKR)
kc_p = [1000, 1000, 1000,1000]; % Pa
% scaling flag
SF = 1;

% % override runout index (number R)
bl_position_nump_or = [1360,1360,1360,1360,1360];
bl_time_or = [3.16,3.60,3.60,3.60,3.60];

% % manually collected runout distances (corresponding to associated work 3)
RA_nump = [nan, nan, 1422, 1450, 1400];
RT_nump = [nan, nan, nan, 1304, 960];
RT_t = [nan, nan, nan, 0.41, 1.26];
RT_dmd0 = [nan, nan,nan, 1.18, 2.18];
RD_nump = [364, 320,364,392];
RD_t = 1E-3*[0.71, 6.75, 66.2, 626];
RD_dmd0 = [496, 107.6, 22.75, 5.42];


%% boulders - impactor mass - decay - rho 3000 - E 50 GPa
% %
%  % file location
%  matfile_arr =[strcat(matfile_location, "decay_boulder_impactorMass/2km_channel_200mfill_mediium_1kms_decay_boulder_impmass0.5_03242024.mat");
%         strcat(matfile_location, "decay_boulder_impactorMass/2km_channel_200mfill_mediium_1kms_decay_boulder_impmass0.75_03072024.mat");
%         strcat(matfile_location, "size_sweep_medium_rho3000_E50GPa_decay/2km_channel_200mfill_medium_1kms_size_1.25m_rho3000_E50GPa_03062024.mat");
%         strcat(matfile_location, "decay_boulder_impactorMass/2km_channel_200mfill_medium_1kms_decay_boulder_impmass1.25_03072024.mat");
%         strcat(matfile_location, "decay_boulder_impactorMass/2km_channel_200mfill_medium_1kms_decay_boulder_impmass2.5_03242024.mat");
%         strcat(matfile_location, "decay_boulder_impactorMass/2km_channel_200mfill_medium_1kms_decay_boulder_impmass5_03242024.mat");];
% 
% Lf = length(matfile_arr);
% evalvals = ["\rho_p/\rho = 0.5", "\rho_p/\rho = 0.75", "\rho_p/\rho= 1", ...
%     "\rho_p/\rho = 1.25", "\rho_p/\rho = 2.5", "\rho_p/\rho = 5"]; % "R = 0.5 mm", "R = 5 mm",
% vel = [1000,1000,1000,1000,1000,1000];
% comptype = 'decay_boulder_impactorMass';
% velyesno = 0;
% r_start_ind = 61;
% r_end_buff = 71;
% tot_deltas_onoff = 0;
% 
% 
% E_p = [50E9,50E9,50E9,50E9,50E9,50E9]; % Pa
% r_p = [1.25,1.25,1.25 ,1.25,1.25,1.25]; % m 0.5E-3, 5E-3
% rho_p = [3000,3000, 3000, 3000, 3000,3000]; % kg/m3
% kc_p = [1000, 1000, 1000, 1000, 1000,1000]; % Pa
% SF = 1;
% 
% % % override runout index (number R)
% bl_position_nump_or = [nan,1360,1360,1360,nan,nan];
% bl_time_or = [nan,3.16,3.16,3.16,nan,nan];
% 
% RA_nump = [nan,1360,1360, 1360,nan,nan];
% RT_nump = [nan,nan, nan, nan,nan,nan];
% RT_t = [nan,nan, nan, nan,nan,nan];
% RT_dmd0 = [nan, nan, nan,nan,nan,nan];
% RD_nump = [267,337,391,433,582,738];
% RD_t = 1E-3*[430,546,626,713,907,1050];
% RD_dmd0 = [5.5,5.5,5.5,5.5,5.5,5.5];


%% boulders - velocity - decay - rho 3000 - E 50 GPa

% % file locatoin
%  matfile_arr =[strcat(matfile_location, "decay_boulder_velocity/2km_channel_200mfill_medium_100ms_boulder_velocity_03072024.mat");
%       strcat(matfile_location, "decay_boulder_velocity/2km_channel_200mfill_medium_177ms_decay_boulder_03042024.mat");
%        strcat(matfile_location, "decay_boulder_velocity/2km_channel_200mfill_medium_316ms_decay_boulder_03042024.mat");
%         strcat(matfile_location, "decay_boulder_velocity/2km_channel_200mfill_medium_562ms_decay_boulder_03042024.mat");
%         strcat(matfile_location, "size_sweep_medium_rho3000_E50GPa_decay/2km_channel_200mfill_medium_1kms_size_1.25m_rho3000_E50GPa_03062024.mat")];
% 
% Lf = length(matfile_arr);
% evalvals = ["v_p = 100 m/s", "v_p = 178 m/s", "v_p = 316 m/s", "v_p = 562 m/s", "v_p = 1 km/s"]; % "R = 0.5 mm", "R = 5 mm",
% vel = [100,178,316,562,1000];
% comptype = 'decay_boulder_velocity';
% velyesno = 0;
% r_start_ind = 61;
% r_end_buff = 71;
% tot_deltas_onoff = 0;
% 
% 
% E_p = [50E9,50E9,50E9,50E9,50E9]; % Pa
% r_p = [1.25,1.25,1.25 ,1.25,1.25]; % m 0.5E-3, 5E-3
% rho_p = [3000,3000, 3000, 3000, 3000]; % kg/m3
% kc_p = [1000, 1000, 1000, 1000, 1000]; % Pa
% SF = 1;
% 
% % % override runout index (number R)
% bl_position_nump_or = [896,1106,1246,1320,1360];
% bl_time_or = [3.16,3.16,3.16,3.16,3.16];
% 
% RA_nump = [135,328,550,952,1150];
% RT_nump = [111,238,450,600,760];
% RT_t = 1E-3*[230,690,976,1209,1619];
% RT_dmd0 = [1.16,1.21,1.12,1.24,1.24];
% RD_nump = [69,126,196,285,391];
% RD_t = 1E-3*[129,218,314,480,626];
% RD_dmd0 = [5.5,5.3,5.6,5.5,5.5];


%% storage and figures
% initialize figures
%h_init_col = figure; axes('NextPlot','add');
% h_init_forcevsradial = figure; axes('NextPlot','add');
%h_init_velocityvsradial = figure; axes('NextPlot','add');
%h_init_VFvsradial = figure; axes('NextPlot','add');
%h_init_row = figure; axes('NextPlot','add');
%h_init_forcevsdepth = figure; axes('NextPlot','add');
%h_init_velocityvsdepth = figure; axes('NextPlot','add');
h_init_VFvsdepth = figure; axes('NextPlot','add');
% h_dil_full = figure; axes('NextPlot','add');
%h_dil_vsdepth = figure; axes('NextPlot','add');
%h_dil_vsradial = figure; axes('NextPlot','add');
% h_dil_partial = figure; axes('NextPlot','add');
h_max_force = figure; axes('NextPlot','add');
%h_max_force_logref = figure; axes('NextPlot','add');
%h_max_velocity = figure; axes('NextPlot','add');
%h_shockdetails_force = figure; axes('NextPlot','add');
%h_shockdetails_runout = figure; axes('NextPlot','add');
% h_shockdetails_speed = figure; axes('NextPlot','add');
% h_shockdetails_speed_exp = figure; axes('NextPlot','add');
% h_shockdetails_speed_arrival = figure; axes('NextPlot','add');
% h_shockdetails_speed_maxforce = figure; axes('NextPlot','add');
%h_dilationbins = figure; axes('NextPlot','add');
h_bulkdilationbinvsdepth = figure; axes('NextPlot','add');
%h_bulkdilationvsradial = figure; axes('NextPlot','add');
%h_heightdilation = figure; axes('NextPlot','add');
%h_heightdilation_max = figure; axes('NextPlot','add');
%h_ss_variances = figure; axes('NextPlot','add');
%h_wave_measurements = figure; axes('NextPlot','add');
%h_init_overlap_v_depth = figure; axes('NextPlot','add');
%h_init_total_overlap_v_depth = figure; axes('NextPlot','add');
%h_sample_height_change_radial = figure; axes('NextPlot','add');
%h_sample_sensor_stress_ratio_yzoyy = figure; axes('NextPlot','add');
%h_density_vs_depth = figure; axes('NextPlot','add');
h_max_overlap = figure; axes('NextPlot','add');
h_max_overlap_vstime = figure; axes('NextPlot','add');
h_max_velocity_new = figure; axes('NextPlot','add');
h_max_overlap_vsdepth = figure; axes('NextPlot','add');
h_max_velocity_new_vsdepth = figure; axes('NextPlot','add');


% saving figures
datetime.setDefaultFormats('default','yyyy-MM-dd')
currenttime = datetime;
formatOut = 'mmddyyyy';
timestring = datestr(currenttime,formatOut);


workingdir = '/Users/hartzell_lab/DOCUMENTS/COLDSPOTS/MATLAB ANALYSIS/ShockInducedDilation/Prod_v/';
plotpath = [workingdir 'ProductionPlots/' comptype '/' timestring '_'];
plotcolor = turbo(Lf);
%plotcolor = flipud(plotcolor);

%harvestpath = [workingdir 'ProcessedData/'];

% storage arrays for steady state dilation error comparison
store_ssdilation = cell(Lf,1);
store_ssdilation_heights = cell(Lf,1);
store_shockrunout = cell(Lf,1);
store_shockspeed = cell(Lf,1);
store_channel = 0;
store_depthplot = 0;
store_initial_avg_phi = cell(Lf,1);
store_initial_std_phi = cell(Lf,1);

store_heights = cell(Lf,1);
store_depths = cell(Lf,1);
store_channels = cell(Lf,1);

store_normalized_channel = cell(Lf,1);
store_dm_normalized = cell(Lf,1);
store_dm_time = cell(Lf,1);

% wave speed -peak forces
store_avg_shockspeed = zeros(Lf,1);
store_std_shockspeed = zeros(Lf,1);

store_avg_shockspeed_exp = zeros(Lf,1);
store_std_shockspeed_exp = zeros(Lf,1);
% wave speed - first arrival
store_avg_shockspeed_arr = zeros(Lf,1);
store_std_shockspeed_arr = zeros(Lf,1);
% wave speed - maxforce
store_avg_shockspeed_max = zeros(Lf,1);
store_std_shockspeed_max = zeros(Lf,1);


% max force
store_max_force = zeros(Lf,1);

store_5_shockspeed = zeros(Lf,1);
store_10_shockspeed = zeros(Lf,1);
store_15_shockspeed = zeros(Lf,1);

ssL = 40; % shock speed length, used for wave speed tracking, vestigial
%store_shockspeeds = zeros(ssL,Lf);
store_shockspeeds = cell(ssL,1);
store_shockspeeds_exp = cell(ssL,1);

store_height_radials = cell(Lf,1);


store_sensor_force_peak_time = zeros(Lf,1);
store_sensor_force_peak = zeros(Lf,1);


% store solitary wave values
% using Pal equation
store_Fm_pal = zeros(Lf,1);
% using our equation (associated work 2)
store_Fm_avg = zeros(Lf,1);
store_Fm_std = zeros(Lf,1);
% using my equation, scaled
store_Fm_scaled_avg = zeros(Lf,1);
store_Fm_scaled_std = zeros(Lf,1);
% simulation value
store_F_SW_avg = zeros(Lf,1);
store_F_SW_std = zeros(Lf,1);
store_F_SW_Fc_avg = zeros(Lf,1);
store_F_SW_Fc_std = zeros(Lf,1);
store_F_SW_Fc_analytic_avg = zeros(Lf,1);
store_F_SW_Fc_analytic_std = zeros(Lf,1);

% store blast loading
store_BL_m = zeros(Lf,1);
store_BL_time = zeros(Lf,1);
store_BL_particleV = zeros(Lf,1);
store_BL_nump = zeros(Lf,1);
store_BL_decayspeed  = zeros(Lf,1);

% sample depth and sample sensor
% pick one depth and report common properties across tests for this depth
% vs channel position (max force vs channel, ss dilation over channel)
% in this depth, pick one radial and for this sensor provide common vs time
% properties (force, stress, dilation, volume fraction)
sample_depth = 4;
sample_radial = 25;
%h_sample_sensor_force = figure; axes('NextPlot','add');
%h_sample_sensor_ystress = figure; axes('NextPlot','add');
%h_sample_sensor_yz_stress = figure; axes('NextPlot','add');
store_peakstress_time = zeros(Lf,1);
store_peakstress_val = zeros(Lf,1);
store_arrival_ind = zeros(Lf,1);
store_arrival_val = zeros(Lf,1);
store_peak_time = zeros(Lf,1);
store_peak_val = zeros(Lf,1);
store_loft_depth = zeros(Lf,1);
% store_peak_time_stress = (Lf,1);
% store_peak_val_stress = (Lf,1);
%h_sample_sensor_stress = figure; axes('NextPlot','add');
%h_sample_sensor_phi = figure; axes('NextPlot','add');
%h_sample_sensor_pdiff_phi = figure; axes('NextPlot','add');
% h_sample_depth_wavemeasurements = figure; axes('NextPlot','add');
h_max_force_new_vsdepth = figure; axes('NextPlot','add');
%h_max_force_vstime = figure; axes('NextPlot','add');
%h_max_force_new_vsdepth_normalized  = figure; axes('NextPlot','add');
%h_max_force_new_vsdepth_scaled  = figure; axes('NextPlot','add');

for j = 1:Lf %1:Lf %1:Lf %:Lf%1:Lf

    % set path varialbes for the investigated .mat file
    matfile_run = matfile_arr(j); % path to .mat file
    evalval_run = convertStringsToChars(evalvals(j)) % label for this test



    %% run the analysis function
    [grainstring, channelstring, porstring, ...
        channel, vy_max, vy_maxtime, percentofmaxvel, Fmag_max, percentofmaxF, ...
        Fmag_colavg_init, Fmag_colstd_init, vy_colavg_init, vy_colstd_init, vf_colavg_init, vf_colstd_init, ...
        depthplot,Fmag_rowavg_init,Fmag_rowstd_init, vy_rowavg_init, vy_rowstd_init, vf_rowavg_init, vf_rowstd_init, ...
        max_forcedist_radial_fullw,max_force_radial_fullw,max_forcestd_radial_fullw, max_forcetime_radial_fullw, ...
        max_vydist_radial_fullw,max_vy_radial_fullw,max_vystd_radial_fullw, max_vytime_radial_fullw, ...
        shockmax,shockmaxstd, shockmin,shockminstd,shockmintime,shock_runout, shock_speed, ...
        shock_runout_F,shock_speed_F,shockmin_F,shockspeedslope, ...
        timeplot,regionstartdepth,regionstopdepth,Fmag_logref_mean,Fmag_logref_std, ...
        bulk_dilation_vsdepth,bulk_dilation_stddv_vsdepth,F_shocktrack_mean,shockspeedslope_arrival,...
        bulk_dilation_vsradial, bulk_dilation_stddv_vsradial, shockspeedslope_maxforce, ...
        height_pdiff_out, ...
        var_vs_channel_x, var_vs_channel_y,sensor_force,sensor_force_arrival,sensor_force_peak,sensor_force_peak_time, ...
        sensor_vfrac,sensor_stressyy,wave_measurement,wave_time,wave_fit,wave_fit_time, ...
        init_overlap_mean_vs_depth,sensor_stress_peak_time,sensor_stress_peak_val,sensor_stress_yz, ...
        store_max_height_pdiff, height_out, shockspeedslope_exp,init_total_overlap_mean_vs_depth,store_heights_bin,radial_sensor_pos, ...
        sensor_overlap, sensor_overlap_peak, sensor_overlap_peak_time,Fmag_shocktrack_maxoverlap,Fmag_shocktrack_maxoverlap_time, ...
        Fmag_shocktrack_maxvelocity,Fmag_shocktrack_maxvelocity_time, F_overlaptrack_mean, F_overlaptrack_mean_time, F_velocitytrack_mean,vf_pdiff_voro, ...
        vf_rowavg_final, vf_rowstd_final, vf_colavg_final, vf_colstd_final,wspacing, F_overlaptrack_mean_vsdepth, F_velocitytrack_mean_vsdepth,F_forcetrack_mean_vsdepth, ...
        F_mag_avg_time,F_mag_avg,loft_depth_track] ...
        = HarvestAndPlotFun_production_loadonly(matfile_run,[workingdir 'Plots/' comptype '/'],evalval_run,comptype,workingdir,vel(j),tot_deltas_onoff,SF,r_p(j));



    %% plot the initial conditions as subplot for the cases


    % string for legend entries
    % look up exponent
    all_vels = [1e-03, 3.1623e-03, 1e-02, 3.1623e-02, 1e-01, 3.1623e-01,...
        5.6234E-01, 1, 1.7783, 3.1623,5.6234, 10];
    all_pow = [-3, -2.5, -2, -1.5, -1, -0.5,-0.25, 0,0.25, 0.5,0.75, 1];
    v_power_ind = find(all_vels == vel(j));
%     initstring = [comptype ' = ' num2str(evalvals(j))];
%     initstring = ['v_{im}  = ' num2str(evalvals(j)) ' m/s'];
    %initstring = ['v_{p}  = 10^{ ' num2str(all_pow(v_power_ind)) '} m/s'];

    initstring = [convertStringsToChars(evalvals(j))]
    % zeros array for horizontal error bars
    znum = length(depthplot);
    zzz = zeros(1,znum);

    dil_start = r_start_ind;
    dil_end = r_end_buff;

    % column averages
%     set(0, 'CurrentFigure', h_init_forcevsradial)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     errorbar(channel./(r_p(j)*100),Fmag_colavg_init,Fmag_colstd_init,'o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:));
%     hold on
%     grid on
%     xlabel('D/R')
%     set(gca, 'YScale', 'log')
%     ylabel('Force [N] (log)')
%     %title(['Initial force vs radial position - ' comptype])
%     set(gca,'fontsize',14)

    %figure
%     set(0, 'CurrentFigure', h_init_velocityvsradial)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     errorbar(channel./(r_p(j)*100),100*vy_colavg_init,100*vy_colstd_init,'--o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:));
%     hold on
%     grid on
%     xlabel('D/R')
%     ylabel('Velocity [cm/s]')
%     %title(['Initial velocity vs radial position - ' comptype])
%     set(gca,'fontsize',14)

%     set(0, 'CurrentFigure', h_init_VFvsradial)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     %figure
%     errorbar(channel./(r_p(j)*100),100*vf_colavg_init,100*vf_colstd_init,'--o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:));
%     hold on
%     grid on
%     xlabel('D/R')
%     ylabel('Volume Fraction [%]')
%     %title(['Initial volume fraction vs radial position - ' comptype])
%     set(gca,'fontsize',14)

    %sgtitle(['Depth averaged initial condition plots - Varied Channel ' comptype])

    % row averages
%     set(0, 'CurrentFigure', h_init_forcevsdepth)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     errorbar(Fmag_rowavg_init,depthplot/(r_p(j)*100),zzz,zzz,Fmag_rowstd_init,Fmag_rowstd_init,'o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:));
%     hold on
%     grid on
%     ylabel('Depth/R')
%     set(gca,'Ydir','reverse')
%     set(gca, 'XScale', 'log')
%     xlabel('Force [N] (log)')
%     %title(['Depth vs Initial force - ' comptype])
%     set(gca,'fontsize',14)

    %figure
%     set(0, 'CurrentFigure', h_init_velocityvsdepth)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     errorbar(100*vy_rowavg_init,depthplot/(r_p(j)*100),zzz,zzz,100*vy_rowstd_init,100*vy_rowstd_init,'--o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:));
%     hold on
%     grid on
%     ylabel('Depth/R')
%     xlabel('Velocity [cm/s]')
%     set(gca,'Ydir','reverse')
%     %title(['Depth vs Initial velocity - ' comptype])
%     set(gca,'fontsize',14)


    %figure
    set(0, 'CurrentFigure', h_init_VFvsdepth)
    %axes('NextPlot','add')
    %subplot(1,1,1)
    errorbar(100*vf_rowavg_init,depthplot/(r_p(j)*100),zzz,zzz,100*vf_rowstd_init,100*vf_rowstd_init,'o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:));
    hold on
    grid on
    ylabel('Depth/R')
    xlabel('\phi_0 [%]')
    set(gca,'Ydir','reverse')
    %title(['Depth vs Initial volume fraction - ' comptype])
    set(gca,'fontsize',14)

    %keyboard
    %sgtitle(['Radially averaged initial condition plots - Varied Channel ' comptype])


%     set(0, 'CurrentFigure', h_init_overlap_v_depth)
%     %axes('NextPlot','add')
%     %subplot(1,1,1)
%     plot(init_overlap_mean_vs_depth,depthplot/(r_p(j)*100),'o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:));
%     hold on
%     grid on
%     xlabel('Coordination # averaged \delta [m]')
%     ylabel('Depth/R')
%     set(gca,'xscale','log')
%     set(gca,'fontsize',14)
%     set(gca,'YDir','reverse')
%     %title('Initial overlaps')
%     legend('show','location','best')

%     set(0, 'CurrentFigure', h_init_total_overlap_v_depth)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     plot(init_total_overlap_mean_vs_depth,depthplot/(100*r_p(j)),'o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:));
%     hold on
%     grid on
%     set(gca,'xscale','log')
%     xlabel('Total particle \delta [m]')
%     ylabel('Depth/R')
%     set(gca,'fontsize',14)
%     set(gca,'YDir','reverse')
%     %title('Initial overlaps')
%     legend('show','location','best')


    %% plot the steady state dilations

    % -- dilation windows -- %
%     set(0, 'CurrentFigure', h_dil_vsradial)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     errorbar(ss_dilation_col_avg_yval./100,ss_dilation_col_avg,ss_dilation_col_std,'--o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     set(gca,'xscale','log')
%     %errorbar(ss_compact_col_avg_yval,ss_compact_col_avg,ss_compact_col_std,'--o','LineWidth',2)
%     xlabel('Channel Position [m]')
%     ylabel('Percent Dilation')
%     %title(['Steady state dilation vs radial position - ' comptype])
%     set(gca,'fontsize',14)

%     zdil = length(ss_dilation_row_avg);
%     zd = zeros(1,zdil);
%     set(0, 'CurrentFigure', h_dil_vsdepth)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     errorbar(ss_dilation_row_avg,ss_dilation_row_avg_yval,zd,zd,ss_dilation_row_std,ss_dilation_row_std,'o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     %errorbar(ss_compact_row_avg_yval,ss_compact_row_avg,ss_compact_row_std,'--o','LineWidth',2)
%     ylabel('Depth [cm]')
%     xlabel('Percent Dilation')
%     set(gca,'Ydir','reverse')
%     %title(['Depth vs Steady state dilation - ' comptype])
%     set(gca,'fontsize',14)

    %sgtitle(['Steady State Dilations - Dilation Windows - Varied Channel ' comptype])

    %% calculate and store loft depth

%     loft_depth_track_ROI = loft_depth_track(:,r_start_ind:r_end_buff);
%     loft_ratio = mean(loft_depth_track_ROI.');
%
%     [loft_max_depth, loft_max_depth_ind] = max(loft_ratio);
%
%     % cut down to just values after max
%     loft_depth_ratio_past_max = loft_ratio(loft_max_depth_ind:end);
%     depth_past_max = depthplot(loft_max_depth_ind:end);
%
%
%     last_loft_depth_ind = find(loft_depth_ratio_past_max < 0.5,1,'first') -1;
%     last_loft_depth = depth_past_max(last_loft_depth_ind);
%     %store_loft_depth(j,1) = last_loft_depth;
%
%     % maybe just use max?
%     last_loft_depth = depthplot(loft_max_depth_ind-1);
%     store_loft_depth(j,1) = last_loft_depth;
%

    [loft_max_depth_pdiff, loft_max_pdiff_depth_ind] = min(vf_pdiff_voro(1:(end-5)));
    % cut down to just values after max
    loft_pdiff_past_max = vf_pdiff_voro(loft_max_pdiff_depth_ind:(end-5));
    depth_past_max = depthplot(loft_max_pdiff_depth_ind:(end-5));

    last_loft_depth_ind = find(loft_pdiff_past_max > -0.5,1,'first');
    if isempty(last_loft_depth_ind)
        last_loft_depth = depthplot(end);
    else
        last_loft_depth = depth_past_max(last_loft_depth_ind);
    end
    store_loft_depth(j,1) = last_loft_depth;

%     to look at while determining algorithm
%     figure
%     plot(depthplot,100*loft_ratio)
%     xlabel('depth [cm]')
%     ylabel('fraction of cells lofting')



    % bulkdilations
    zdilb = length(bulk_dilation_vsdepth);
    zdb = zeros(1,zdilb);
    set(0, 'CurrentFigure', h_bulkdilationbinvsdepth)
    %axes('NextPlot','add')
    %subplot(1,1,1)
    errorbar(vf_pdiff_voro(1:(end-2)),depthplot(1:(end-2))/(100*r_p(j)),zdb(1:(end-2))/(100*r_p(j)),zdb(1:(end-2)),bulk_dilation_stddv_vsdepth(1:(end-2)),bulk_dilation_stddv_vsdepth(1:(end-2)),'o:','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:))
    hold on
    grid on

    % find index of loft depth
    [~,lld_ind] = min(abs(depthplot-last_loft_depth));
    if isempty(lld_ind)
        lld_ind = length(depthplot);
    end
    %plot(vf_pdiff_voro(lld_ind),last_loft_depth/(100*r_p(j)),'o','LineWidth',2,'Color',plotcolor(j,:),'HandleVisibility','off','MarkerSize',20)
    %errorbar(ss_compact_row_avg_yval,ss_compact_row_avg,ss_compact_row_std,'--o','LineWidth',2)
    ylabel('Depth/R')
    xlabel('%\Delta\rho_{bulk}')
    %ylim([0, 25])
    set(gca,'Ydir','reverse')
    %title(['Depth vs Bulk Dilation - ' comptype])
    set(gca,'fontsize',14)



    %bulk_dilation_vsradial, bulk_dilation_stddv_vsradial
%     set(0, 'CurrentFigure', h_bulkdilationvsradial)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     semilogx(channel./(r_p(j)*100),bulk_dilation_vsradial,'--o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     %errorbar(ss_compact_row_avg_yval,ss_compact_row_avg,ss_compact_row_std,'--o','LineWidth',2)
%     ylabel('Percent Dilation')
%     xlabel('D/R')
%     %set(gca,'Ydir','reverse')
%     title(['New Bulk Dilation vs Radial - ' comptype])
%     set(gca,'fontsize',14)

%% sample sensor
    % plot sensor values - vfrac
%     set(0, 'CurrentFigure', h_sample_sensor_phi);
%     axes('NextPlot','add');
%     subplot(1,1,1)
%     plot(timeplot,sensor_vfrac,'.--','LineWidth',2','DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     legend('show','location','best')
%     %title(['\phi - time, ' ', depth ind ' num2str(sample_depth) ', radial ind ' num2str(sample_radial) ', ' comptype])
%
%     set(gca,'fontsize',14)
%     xlabel('Time [s]')
%     ylabel('\phi')

    % plot sensor values - %diff vfrac
%     set(0, 'CurrentFigure', h_sample_sensor_pdiff_phi);
%     axes('NextPlot','add');
%     subplot(1,1,1)
%     plot(timeplot,(sensor_vfrac-sensor_vfrac(1))/sensor_vfrac(1),'.','LineWidth',2','DisplayName',[initstring],'Color',plotcolor(j,:))
%     legend('show','location','best')
%     title(['\phi - time, ' ', depth ind ' num2str(sample_depth) ', radial ind ' num2str(sample_radial) ', ' comptype])
%
%     set(gca,'fontsize',14)
%     xlabel('Time [s]')
%     ylabel('% difference \phi')

 %% dilation
%     % dilation as end bed height difference
%     set(0, 'CurrentFigure', h_heightdilation)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     plot(channel/(r_p(j)*100),height_out/r_p(j),'o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     ylabel('Change in height / R ')
%     xlabel('D/R')
%     legend('show','location','best')
%     %title(['Bed Height % difference SS'])
%     %ylim([-1,2])
%     set(gca,'fontsize',14)

%     set(0, 'CurrentFigure', h_ss_variances)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     hold on
%     plot(var_vs_channel_x,var_vs_channel_y,'o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:))
%     %semilogy(channel(jj)/100,abs(mean(vfrac_vs_time_mean_rad(jj,(end-10):end))),'ko')
%     grid on
%     xlabel('R [m]')
%     ylabel('10^3\sigma_{\phi}')
%     %title(['R < 30 cm, V_{im} = ' num2str(evalval), ' [m/s]'])
%     title(['Average \sigma\phi_{ss}, ' comptype])
%     set(gca,'fontsize',14)
%     set(gca, 'YScale', 'log')

%     set(0, 'CurrentFigure', h_heightdilation_max)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     plot(channel/(r_p(j)*100),store_max_height_pdiff,'o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     ylabel('% difference z')
%     xlabel('D/R')
%     legend('show','location','best')
%     title(['Bed Height % difference MAX'])
%     %ylim([-2,20])


    % new figure that plots our densitites vs Bandfield densities
    rho_p_gcc = 2.5; % enter particle bulk density that was input into sim g/cc
%     set(0, 'CurrentFigure', h_density_vs_depth)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     vf_rowavg_init(end-1) = NaN;
%     vf_rowavg_final(end-1) = NaN;
%     plot(rho_p_gcc*vf_rowavg_init,depthplot/(100*r_p(j)),'o:','LineWidth',2,'DisplayName',['\rho_{0, ' initstring '}'],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     plot(rho_p_gcc*vf_rowavg_final,depthplot/(100*r_p(j)),'+-','LineWidth',2,'DisplayName',['\rho_{f, ' initstring '}'],'Color',plotcolor(j,:))
%     xlabel('Bulk density [g/cm^3]')
%     ylabel('Depth/R')
%     %ylim([0, 25])


    %% max velocity and force (averaged
    % channel length (force radial)
%     set(0, 'CurrentFigure', h_max_force_logref)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     % semilogy(max_forcedist_radial_dw,max_force_radial_dw,'--o','LineWidth',2)
%     %errorbar(max_forcedist_radial_fullw,max_force_radial_fullw,max_forcestd_radial_fullw,':o','LineWidth',2,'DisplayName',[initstring])
% %     errorbar(channel,Fmag_logref_mean,Fmag_logref_std,':o','LineWidth',2,'DisplayName',[initstring])
%     plot(channel/(r_p(j)*100),Fmag_logref_mean,'o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:))
%     grid on
%     hold on
%     set(gca, 'YScale', 'log')
%     xlabel('D/R')
%     ylabel('Force [N]')
%     %hold off
%     %yyaxis right
%     %plot(max_forcedist_radial_fullw,max_forcetime_radial_fullw,'DisplayName',[initstring])
%     %ylabel('Time [s]')
%     title(['Averaged force compared to the initial force OOM vs channel position - ' comptype])
    %legend('Max force', 'Time of max force','Location','Best')



    set(0, 'CurrentFigure', h_max_force)
    axes('NextPlot','add')
    subplot(1,1,1)
    % semilogy(max_forcedist_radial_dw,max_force_radial_dw,'--o','LineWidth',2)
    %errorbar(max_forcedist_radial_fullw,max_force_radial_fullw,max_forcestd_radial_fullw,':o','LineWidth',2,'DisplayName',[initstring])
%     errorbar(channel,Fmag_logref_mean,Fmag_logref_std,':o','LineWidth',2,'DisplayName',[initstring])
    plot(channel/(r_p(j)*100),F_shocktrack_mean,'o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:))
    grid on
    hold on
    set(gca, 'YScale', 'log')
    xlabel('D/R')
    ylabel('Force [N]')
    %hold off
    %yyaxis right
    %plot(max_forcedist_radial_fullw,max_forcetime_radial_fullw,'DisplayName',[initstring])
    %ylabel('Time [s]')
    %title(['Averaged force vs channel position - ' comptype])
    %legend('Max force', 'Time of max force','Location','Best')
%     % ** velocity **

    % channel length (velocity radial)
%     set(0, 'CurrentFigure', h_max_velocity)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     errorbar(max_vydist_radial_fullw/(r_p(j)*100),100*max_vy_radial_fullw,100*max_vystd_radial_fullw,'--o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     %plot(max_vydist_radial_fullw,100*max_vy_radial_fullw,':o','LineWidth',2)
%     xlabel('D/R')
%     ylabel('Velocity [cm/s]')
%     %hold off
%     %yyaxis right
%     %plot(max_vydist_radial_fullw, max_vytime_radial_fullw,'DisplayName',[initstring])
%     %title(['Averaged maximum velocity vs channel position - ' comptype])
%     %legend('Max velocity','Time of max velocity', 'All windows','Location','Best')
%
%     %sgtitle(['Maximum Force and Velocity vs Radial position for varied ' comptype])
%
    %% find the length of the blast loading region
    % density multiplier - accounts for poured/impacted channel density
    % difference
%     rho_multi = 2500/rho_p(j);
    % find the peak
    l_id = length(init_overlap_mean_vs_depth); % to match region used in finding max dm
    dmd0 = F_overlaptrack_mean./nanmean(init_overlap_mean_vs_depth(regionstartdepth:regionstopdepth));
%     if r_p(j) == 12.5E-3
%         dmd0 = dmd0*rho_multi;
%     end
    [dmd0_max,dmd0max_ind] = max(dmd0);

    % find the average from the analysis region
    dmd0_avg = nanmean(dmd0(dil_start:dil_end));
    dmd0_std = nanstd(dmd0(dil_start:dil_end));

    % cut arrays down to just be max forward
    dmd0_forward = dmd0(dmd0max_ind:end);
    channel_forward = channel(dmd0max_ind:end);
    time_forward = F_overlaptrack_mean_time(dmd0max_ind:end);
    vel_forward = F_velocitytrack_mean(dmd0max_ind:end);

    % find the first time the overlap returns to be within 1 std dev of avg
    % blast loading index
    %bl_ind = find(dmd0_forward < dmd0_max*0.20,1,'first');
    %bl_ind = find(dmd0_forward < (dmd0_avg+3*dmd0_std),1,'first');
    bl_ind = find(dmd0 < 1.01,1,'first');

    if isempty(bl_ind)
        %keyboard
        bl_position_meters = nan;
        bl_position_nump = nan;
        bl_position_time = nan;
        bl_position_particle_vel = nan;
    else
        bl_position_meters = channel_forward(bl_ind)/100;
        bl_position_nump = (bl_position_meters)/(r_p(j));
        bl_position_time = time_forward(bl_ind);
        bl_position_particle_vel = vel_forward(bl_ind);
    end
    store_BL_m(j,1) = bl_position_meters;
    store_BL_nump(j,1) = bl_position_nump;
    store_BL_time(j,1) = bl_position_time;
    store_BL_particleV(j,1) = bl_position_particle_vel;
    store_BL_decayspeed(j,1) = bl_position_meters/bl_position_time;

    set(0, 'CurrentFigure', h_max_overlap)
    %axes('NextPlot','add')
    %subplot(1,1,1)
    plot(channel/(r_p(j)*100),dmd0,'o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:))
    grid on
    hold on
    %if j > 1
%     plot(bl_position_nump,dmd0(bl_ind),'ko','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')
% %     plot(bl_position_nump,dmd0(bl_ind),'kx','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')
%         plot(bl_position_nump_or(j),1.05,'ko','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')
%     plot(bl_position_nump_or(j),1.05,'kx','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')
    % RA
    %plot(RA_nump(j),1.05,'ro','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')
    %plot(RA_nump(j),1.05,'rx','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')
    % RT
    %plot(RT_nump(j),RT_dmd0(j),'r+','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')
    %plot(RT_nump(j),RT_dmd0(j),'ro','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')
    % RD
    plot(RD_nump(j),RD_dmd0(j),'ro','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')
%     plot(RD_nump(j),RD_dmd0(j),'','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    xlabel('D/R')
    ylabel('\delta_m/\delta_0')
    set(gca,'fontsize',16)
    legend('show','location','SouthWest')





    % fit and plot
    % varied particle sizes
    xvar_channel = channel/(r_p(j)*100); % non dimensional


    % need to omit first point (x=0) for power law fit
    xvar_channel = xvar_channel(2:end);
    dmd0_fit = dmd0(2:end);

    myfun = @(x,xdata) x(1)*xdata.^(x(2)) + 1;
    %myfun = @(x,xdata) x(1)*exp(x(2)*xdata);
    p0 = [30,-0.7055]; % initial guess

   % [p,resnorm,residual,exitflag,output,lambda,jacobian] =  ...
   %     lsqcurvefit(myfun,p0,xvar_channel,dmd0_fit);

    % each row gives the upper and lower 95% confidence bounds
  %  conf = nlparci(p,residual,'jacobian',jacobian);
    %p_up = conf(1:2,1);
    %p_low = conf(1:2,2);


    %size_uniform = logspace(-5,2,50);
   % gen_uniform = myfun(p,xvar_channel);

    %plot(xvar_channel,gen_uniform,'g-','LineWidth',2,'HandleVisibility','off')



    %%
    set(0, 'CurrentFigure', h_max_overlap_vstime)
    %axes('NextPlot','add')
    %subplot(1,1,1)
    plot(F_overlaptrack_mean_time,dmd0,'o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:))
    grid on
    hold on
    %if j > 1
    %plot(bl_position_nump,dmd0(bl_ind),'ko','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')
    %plot(bl_position_nump,dmd0(bl_ind),'kx','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')
    %end
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    xlabel('Time (s)')
    ylabel('\delta_m/\delta_0')
    set(gca,'fontsize',16)
    legend('show','location','SouthWest')
%     plot(bl_time_or(j),1.05,'ko','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')
%     plot(bl_time_or(j),1.05,'kx','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')
%

  % RT
    %plot(RT_t(j),RT_dmd0(j),'r+','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')
    %plot(RT_t(j),RT_dmd0(j),'ro','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')
    % RD
    plot(RD_t(j),RD_dmd0(j),'ro','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')


    scaler = (1/1)*5*mean(F_overlaptrack_mean_time.^(-0.71));
    %plot(F_overlaptrack_mean_time,scaler*F_overlaptrack_mean_time.^(-0.71),'k--','LineWidth',2,'HandleVisibility','off')


    % store all the overlap things here
    store_normalized_channel{j,1} = channel/(r_p(j)*100);
    store_dm_normalized{j,1} = dmd0;
    store_dm_time{j,1} = F_overlaptrack_mean_time;

    %%
    set(0, 'CurrentFigure', h_max_velocity_new)
    %axes('NextPlot','add')
    %subplot(1,1,1)
    % semilogy(max_forcedist_radial_dw,max_force_radial_dw,'--o','LineWidth',2)
    %errorbar(max_forcedist_radial_fullw,max_force_radial_fullw,max_forcestd_radial_fullw,':o','LineWidth',2,'DisplayName',[initstring])
%     errorbar(channel,Fmag_logref_mean,Fmag_logref_std,':o','LineWidth',2,'DisplayName',[initstring])
    plot(channel/(r_p(j)*100),F_velocitytrack_mean,'o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:))

    c_nd = channel/(r_p(j)*100);
    [~,velt_ind] = min(abs(c_nd-RT_nump(j)));

    % RT
   % plot(RT_nump(j),F_velocitytrack_mean(velt_ind),'r+','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')
   % plot(RT_nump(j),F_velocitytrack_mean(velt_ind),'ro','LineWidth',2,'MarkerSize',25,'HandleVisibility','off')

    grid on
    hold on
    yline(4,'LineWidth',2,'HandleVisibility','off')
    set(gca, 'YScale', 'log')
    set(gca, 'XScale', 'log')
    xlabel('D/R')
    ylabel('v_m [m/s]')
    set(gca,'fontsize',16)
    legend('show','Location','SouthWest')






    set(0, 'CurrentFigure', h_max_overlap_vsdepth)
    %axes('NextPlot','add')
    %subplot(1,1,1)
    plot(F_overlaptrack_mean_vsdepth,depthplot/(100*r_p(j)),'+','LineWidth',2','DisplayName',['\delta_m, ' initstring],'Color',plotcolor(j,:))
    hold on
    grid on
    if 1% j == Lf
        plot(init_overlap_mean_vs_depth,depthplot/(100*r_p(j)),'o','LineWidth',2,'DisplayName',['\delta_0, ' initstring],'Color',plotcolor(j,:),'HandleVisibility','off');
    end
    set(gca,'Ydir','reverse')
    ylabel('Depth/R')
    xlabel('\delta [m]')
    %title(['Depth vs Wave speed - ' comptype])
    %title('peak force')
    legend('show','Location','best')
    set(gca,'fontsize',14)
    set(gca,'xscale','log')


    % store avereage initial phi
    % omit first row and bottom 1cm%
    store_initial_avg_phi{j,1} = 100*nanmean(vf_rowavg_init((regionstartdepth+1):regionstopdepth));
    store_initial_std_phi{j,1} = 100*nanstd(vf_rowstd_init((regionstartdepth+1):regionstopdepth));

    %keyboard
    [init_overlap_calc,init_overlap_calc_cohesion] = ...
        init_overlap_calc_fun(rho_p(j),cell2mat(store_initial_avg_phi(j))/100, ...
        depthplot/100,r_p(j),1.625,E_p(j),kc_p(j),0.2,vel(j));
    % calculate initial overlap
    if 1% j == Lf
        plot(init_overlap_calc,depthplot/(100*r_p(j)),'--','LineWidth',2,'HandleVisibility','off','Color',plotcolor(j,:));
    end


    set(0, 'CurrentFigure', h_max_velocity_new_vsdepth)
    %axes('NextPlot','add')
    %subplot(1,1,1)
    plot(F_velocitytrack_mean_vsdepth,depthplot/(100*r_p(j)),'o','LineWidth',2','DisplayName',[initstring],'Color',plotcolor(j,:))
    hold on
    grid on
    set(gca,'Ydir','reverse')
    ylabel('Depth/R')
    xlabel('v_m [m/s]')
    %title(['Depth vs Wave speed - ' comptype])
    %title('peak force')
    legend('show','Location','best')
    set(gca,'fontsize',14)

  set(0, 'CurrentFigure', h_max_force_new_vsdepth)
    %axes('NextPlot','add')
    %subplot(1,1,1)
    plot(F_forcetrack_mean_vsdepth,depthplot/(100*r_p(j)),'o','LineWidth',2','DisplayName',[initstring],'Color',plotcolor(j,:))
    hold on
    grid on
    set(gca,'Ydir','reverse')
    ylabel('Depth/R')
    xlabel('F_m [N]')
    %title(['Depth vs Wave speed - ' comptype])
    %title('peak force')
    legend('show','Location','best')
    set(gca,'fontsize',14)
    set(gca,'xscale','log')

    pvelocity_SW = mean(F_velocitytrack_mean(dil_start:dil_end));
    pvelocity_SW_std = std(F_velocitytrack_mean(dil_start:dil_end));

%     [Fm] = FSW_impulse_fun(F_velocitytrack_mean_vsdepth.',init_overlap_mean_vs_depth);
%     [Fm_vs_depth] = FSW_impulse_fun(10,init_overlap_mean_vs_depth);

%     r_p = 1.25E-3;
%     E_p = 5E6;
%     rho_p = 2500;


%     cell2mat(store_initial_avg_phi)
%     [Fm_vs_depth,Fm_vs_depth_scaled, ...
%         Fm_pal, Fm_avg, Fm_std, Fm_scaled_avg, Fm_scaled_std,Fc,Fc_analytical] ...
%         = FSW_impulse_fun_v4(vel(j),init_overlap_mean_vs_depth,1/100, ...
%         E_p(j),r_p(j),rho_p(j),depthplot,store_initial_avg_phi{j,1});

    % compute the predicted force in the wavefront using the equation from
    % associated work 2
     [Fm_vs_depth,Fm_vs_depth_scaled, ...
        Fm_pal, Fm_avg, Fm_std, Fm_scaled_avg, Fm_scaled_std,Fc,Fc_analytical] ...
        = FSW_impulse_fun_v4(vel(j),init_overlap_calc_cohesion,1/100, ...
        E_p(j),r_p(j),rho_p(j),depthplot,store_initial_avg_phi{j,1},kc_p(j));


    %plot(Fm_vs_depth,depthplot/(100*r_p(j)),'*','LineWidth',2','DisplayName',[initstring ': predicted'],'Color',plotcolor(j,:))
    plot(Fm_vs_depth,depthplot/(100*r_p(j)),'+','LineWidth',2','DisplayName',[initstring ': predicted'],'Color',plotcolor(j,:),'HandleVisibility','off')

    % using Pal equation
    store_Fm_pal(j,1) = Fm_pal;
    % using my equation
    store_Fm_avg(j,1) = Fm_avg;
    store_Fm_std(j,1) = Fm_std;
    % using my equation, scaled
    store_Fm_scaled_avg(j,1) = Fm_scaled_avg;
    store_Fm_scaled_std(j,1) = Fm_scaled_std;
    % simulation value - maybe remove top and bottom -most depths from calc?
    store_F_SW_avg(j,1) = mean(F_forcetrack_mean_vsdepth);
    store_F_SW_std(j,1) = std(F_forcetrack_mean_vsdepth);
    store_F_SW_Fc_avg(j,1) = mean(Fc);
    store_F_SW_Fc_std(j,1) = std(Fc);
    store_F_SW_Fc_analytic_avg(j,1) = mean(Fc_analytical);
    store_F_SW_Fc_analytic_std(j,1) = std(Fc_analytical);

    %keyboard

   % keyboard

%     set(0, 'CurrentFigure', h_max_force_new_vsdepth_normalized)
%     %axes('NextPlot','add')
%     %subplot(1,1,1)
%     %plot(F_forcetrack_mean_vsdepth/F_forcetrack_mean_vsdepth(length(depthplot)-3),depthplot/(100*r_p(j)),'o','LineWidth',2','DisplayName',[initstring ': sim'],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     %plot(Fm_vs_depth/Fm_vs_depth(length(depthplot)-3),depthplot/(100*r_p(j)),'*','LineWidth',2','DisplayName',[initstring ': predicted'],'Color',plotcolor(j,:))
%     set(gca,'Ydir','reverse')
%     ylabel('Depth/R')
%     xlabel('F_m/F_m(floor)')
%     %title(['Depth vs Wave speed - ' comptype])
%     %title('peak force')
%     legend('show','Location','best')
%     set(gca,'fontsize',14)
%     set(gca,'xscale','log')

%     set(0, 'CurrentFigure', h_max_force_new_vsdepth_scaled)
%     %axes('NextPlot','add')
%     %subplot(1,1,1)
%     plot((F_forcetrack_mean_vsdepth.')./Fm_vs_depth,depthplot/(100*r_p(j)),'o','LineWidth',2','DisplayName',[initstring ': sim/pred'],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     %plot(Fm/Fm(length(depthplot)-3),depthplot,'*','LineWidth',2','DisplayName',[initstring ': predicted'],'Color',plotcolor(j,:))
%     set(gca,'Ydir','reverse')
%     ylabel('Depth/R')
%     xlabel('F_{m,sim}/F_{m,pred} [N]')
%     %title(['Depth vs Wave speed - ' comptype])
%     %title('peak force')
%     legend('show','Location','best')
%     set(gca,'fontsize',14)
%     set(gca,'xscale','log')


    %% Plot force shock details
%     zshock = length(regionstartdepth:regionstopdepth);
%     zs = zeros(1,zshock);
%     set(0, 'CurrentFigure', h_shockdetails_force)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     errorbar(shockmax(regionstartdepth:regionstopdepth),depthplot(regionstartdepth:regionstopdepth)/(100*r_p(j)),zs,zs,shockmaxstd(regionstartdepth:regionstopdepth),shockmaxstd(regionstartdepth:regionstopdepth),'--o','LineWidth',2,'DisplayName',['Max - ' channelstring 'm c, ' grainstring ' cm f'])
%     hold on
%     grid on
%     errorbar(shockmin_F(regionstartdepth:regionstopdepth),depthplot(regionstartdepth:regionstopdepth)/(100*r_p(j)),zs,zs,shockminstd(regionstartdepth:regionstopdepth),shockminstd(regionstartdepth:regionstopdepth),'--x','LineWidth',2,'DisplayName',['Min - ' channelstring 'm c, ' grainstring ' cm f'])
%     set(gca, 'XScale', 'log')
%     set(gca,'Ydir','reverse')
%     xlabel('Force [log(N)]')
%     ylabel('Depth/R')
%     legend('show','Location','best')
%     %title(['Depth vs Wave maximum and minimum forces - ' comptype])
%     %legend('Max','Min')
%     set(gca,'fontsize',14)


%     set(0, 'CurrentFigure', h_shockdetails_runout)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     plot(shock_runout_F,depthplot/(100*r_p(j)),'--o','LineWidth',2,'DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     ylabel('Depth/R')
%     xlabel('Radial distance [m]')
%     set(gca,'Ydir','reverse')
%     %title(['Depth vs Wave runout distance - ' comptype])
%     legend('show','Location','best')
%     set(gca,'fontsize',14)



%     set(0, 'CurrentFigure', h_shockdetails_speed)
%     %axes('NextPlot','add')
%     %subplot(1,1,1)
%     plot(shockspeedslope,depthplot/(100*r_p(j)),'o','LineWidth',2','DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     set(gca,'Ydir','reverse')
%     ylabel('Depth/R')
%     xlabel('c_w [m/s]')
%     %title(['Depth vs Wave speed - ' comptype])
%     %title('peak force')
%     legend('show','Location','best')
%     set(gca,'fontsize',14)


%     set(0, 'CurrentFigure', h_shockdetails_speed_exp)
%     %axes('NextPlot','add')
%     %subplot(1,1,1)
%     plot(shockspeedslope_exp,depthplot/(100*r_p(j)),'o','LineWidth',2','DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     set(gca,'Ydir','reverse')
%     ylabel('Depth/R')
%     xlabel('c_w [m/s]')
%     %title(['Depth vs Wave speed - ' comptype])
%     title('peak force - exp')
%     legend('show','Location','best')
%     set(gca,'fontsize',14)

%     set(0, 'CurrentFigure', h_shockdetails_speed_arrival)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     plot(shockspeedslope_arrival,depthplot/(100*r_p(j)),'--o','LineWidth',2','DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     set(gca','Ydir','reverse')
%     ylabel('Depth/R')
%     xlabel('Wave Speed [m/s]')
%     %title(['Depth vs Wave speed - ' comptype])
%     %title('arrival')
%     legend('show','Location','best')
%     set(gca,'fontsize',14)
%     %sgtitle(['Shock strength, duration, speed and runout vs depth using Force for varied channel ' comptype ])


%     set(0, 'CurrentFigure', h_shockdetails_speed_maxforce)
%     axes('NextPlot','add')
%     subplot(1,1,1)
%     plot(shockspeedslope_maxforce,depthplot/(100*r_p(j)),'--o','LineWidth',2','DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     set(gca','Ydir','reverse')
%     ylabel('Depth/R')
%     xlabel('Velocity [m/s]')
%     title('max force')
%     %title(['Depth vs Wave speed - ' comptype])
%     legend('show','Location','best')
%     set(gca,'fontsize',14)
%     %sgtitle(['Shock strength, du
%
 %% sensor plotting


    % plot sensor values - force
%     set(0, 'CurrentFigure', h_sample_sensor_force);
%     axes('NextPlot','add');
%     subplot(1,1,1)
%     semilogy(timeplot,abs(sensor_force),'-','LineWidth',2','DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     %title(['Force v - time, ' ', depth ind ' num2str(sample_depth) ', radial ind ' num2str(sample_radial) ', ' comptype])
%     %semilogy(sensor_force_peak_time,abs(sensor_force_peak),'ro','LineWidth',2','DisplayName',['Peak'])
%     legend('show','location','best')
%     set(gca,'fontsize',14)
%     xlabel('Time [s]')
%     ylabel('Average Particle Force [N]')

%     set(0, 'CurrentFigure', h_sample_sensor_ystress);
%     axes('NextPlot','add');
%     subplot(1,1,1)
%     semilogy(timeplot,abs(sensor_stressyy),'.','LineWidth',2','DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     %title(['Stress yy v time, ' ', depth ind ' num2str(sample_depth) ', radial ind ' num2str(sample_radial) ', ' comptype])
%     legend('show','location','best')
%     set(gca,'fontsize',14)
%     xlabel('Time [s]')
%     ylabel('Average Particle Stress - YY [N*m]')

%     set(0, 'CurrentFigure', h_sample_sensor_yz_stress);
%     axes('NextPlot','add');
%     subplot(1,1,1)
%     semilogy(timeplot,abs(sensor_stress_yz),'.','LineWidth',2','DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     %title(['Stress yy v time, ' ', depth ind ' num2str(sample_depth) ', radial ind ' num2str(sample_radial) ', ' comptype])
%     legend('show','location','best')
%     set(gca,'fontsize',14)
%     xlabel('Time [s]')
%     ylabel('Average Particle Stress - YZ [N*m]')

%     set(0, 'CurrentFigure', h_sample_sensor_stress_ratio_yzoyy);
%     axes('NextPlot','add');
%     subplot(1,1,1)
%     semilogy(timeplot,abs(sensor_stress_yz./sensor_stressyy),'.--','LineWidth',2','DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     %title(['Stress yy v time, ' ', depth ind ' num2str(sample_depth) ', radial ind ' num2str(sample_radial) ', ' comptype])
%     legend('show','location','best')
%     set(gca,'fontsize',14)
%     xlabel('Time [s]')
%     ylabel('\tau_{yz}/\tau_{yy}')

    % plot sensor height change vs time
%     set(0, 'CurrentFigure', h_sample_height_change_radial);
%     axes('NextPlot','add');
%     subplot(1,1,1)
%     semilogy(timeplot,(store_heights_bin - store_heights_bin(1))/(r_p(j)),'o','LineWidth',2','DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     title(['Height change at ' num2str(radial_sensor_pos) ' m'])
%     legend('show','location','best')
%     set(gca,'fontsize',14)
%     xlabel('Time [s]')
%     ylabel('Height/R')


    % plot sensor values - vfrac
%     set(0, 'CurrentFigure', h_sample_sensor_phi);
%     axes('NextPlot','add');
%     subplot(1,1,1)
%     plot(timeplot,sensor_vfrac,'.','LineWidth',2','DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     legend('show','location','best')
%     %title(['\phi - time, ' ', depth ind ' num2str(sample_depth) ', radial ind ' num2str(sample_radial) ', ' comptype])
%
%     set(gca,'fontsize',14)
%     xlabel('Time [s]')
%     ylabel('\phi')

    % plot sensor values - %diff vfrac
%     set(0, 'CurrentFigure', h_sample_sensor_pdiff_phi);
%     axes('NextPlot','add');
%     subplot(1,1,1)
%     plot(timeplot,(sensor_vfrac-sensor_vfrac(1))/sensor_vfrac(1),'.','LineWidth',2','DisplayName',[initstring],'Color',plotcolor(j,:))
%     legend('show','location','best')
%     title(['\phi - time, ' ', depth ind ' num2str(sample_depth) ', radial ind ' num2str(sample_radial) ', ' comptype])

%     set(gca,'fontsize',14)
%     xlabel('Time [s]')
%     ylabel('% difference \phi')

    % plot sensor values - %wave measurements
%     set(0, 'CurrentFigure', h_wave_measurements);
%     axes('NextPlot','add');
%     subplot(1,1,1)
%     plot(wave_time,wave_measurement/r_p(j),'o','LineWidth',2','DisplayName',[initstring],'Color',plotcolor(j,:))
%     hold on
%     grid on
%     %plot(wave_fit_time,wave_fit,'--','LineWidth',2','DisplayName',['Fit, ' initstring],'Color',plotcolor(j,:))
%     set(gca,'fontsize',14)
%     xlabel('Time [s]')
%     ylabel('Sensor location / R')
%     legend('show','location','best')

    %keyboard

    % store the longest channel length or depthplot
    if length(channel) > length(store_channel)
        store_channel = channel;
    end
    if length(depthplot) > length(store_depthplot)
        store_depthplot = depthplot;
    end

    % store arrival and peaks
    store_arrival_ind(j,1) = sensor_force_arrival;
    store_arrival_val(j,1) = sensor_force(sensor_force_arrival);
    store_peak_time(j,1) = sensor_force_peak_time;
    store_peak_val(j,1) = sensor_force_peak;
    store_peakstress_time(j,1) = sensor_stress_peak_time;
    store_peakstress_val(j,1) = sensor_stress_peak_val;
%     store_peak_time_stress(j,1) = ;
%     store_peak_val_stress(j,1) = ;

    % store dilations
    %store_ssdilation{j,1} = dilationpdiff_SS;
    store_ssdilation_heights{j,1} = height_pdiff_out;
    % define range over which to find average dilation
    dil_start = r_start_ind;% 1m
    dil_end = r_end_buff; % 1m, 50cm

    %for length test only
%     if j == 4
%         dil_start = 51;  % 1m
%         dil_end = 100; % 1m, 50cm
%     end

    if j == 1
        channel_old = channel;
        r_p_old = r_p(j);
        %keyboard % send back to j == 1 after plot making

    end

%     store_ssdilaiton_hiehgt_mean(j,1) = mean(height_pdiff_out(dil_start:dil_end));
%     store_ssdilaiton_hiehgt_std(j,1) = std(height_pdiff_out(dil_start:dil_end));
    % compute bounds based on window spacing - want to chop off 50 cm in
    % both dir

%     % use below when particle sizes are changing
%     new_dil_start = floor(0.5/wspacing);
%     new_dil_end = length(channel) - new_dil_start;
    new_dil_start = dil_start;
    new_dil_end = dil_end;
    % use this to select a specific region in the 2m, 2.5mm diameter case
%     new_dil_start = r_start_ind;
%     new_dil_end = r_end_buff;
%
    %keyboard
    store_ssdilation_height_mean(j,1) = 100*mean(height_out(new_dil_start:new_dil_end));
    store_ssdilation_height_std(j,1) = 100*std(height_out(new_dil_start:new_dil_end));

    store_solitary_overlap_mean(j,1) = mean(F_overlaptrack_mean(dil_start:dil_end));
    store_solitary_overlap_std(j,1) = std(F_overlaptrack_mean(dil_start:dil_end));

    store_solitary_pvelocity_mean(j,1) = mean(F_velocitytrack_mean(dil_start:dil_end));
    store_solitary_pvelocity_std(j,1) = std(F_velocitytrack_mean(dil_start:dil_end));

    store_init_overlap_mean(j,1) = nanmean(init_overlap_mean_vs_depth);
    store_init_overlap_std(j,1) = nanstd(init_overlap_mean_vs_depth);

    % store shock details
    store_shockrunout{j,1} = shock_runout_F;
    store_shockspeed{j,1} = shockspeedslope;
    store_shockspeeds_exp{j,1} = shockspeedslope_exp;
    store_depths{j,1} = depthplot;

    % start a few cm below surface since surface is noisey
    wavestart = 3;
    % stop approx 5 particle diameter from wall to account for wall sheath
    % size (gupta, 2016)
    % 5 particle diamter = 5*2.5mm = 10.5mm = 12.5mm = 1.25 cm = ~3 grids
    waveend = length(shockspeedslope_arrival) - 3;

    %keyboard
    store_avg_shockspeed(j,1) = nanmean(shockspeedslope(wavestart:waveend));
    store_std_shockspeed(j,1) = nanstd(shockspeedslope(wavestart:waveend));

    store_avg_shockspeed_exp(j,1) = nanmean(shockspeedslope_exp(wavestart:waveend));
    store_std_shockspeed_exp(j,1) = nanstd(shockspeedslope_exp(wavestart:waveend));

    store_avg_shockspeed_arr(j,1) = nanmean(shockspeedslope_arrival(wavestart:waveend));
    store_std_shockspeed_arr(j,1) = nanstd(shockspeedslope_arrival(wavestart:waveend));


    store_avg_shockspeed_max(j,1) = nanmean(shockspeedslope_maxforce(wavestart:waveend));
    store_std_shockspeed_max(j,1) = nanstd(shockspeedslope_maxforce(wavestart:waveend));

    store_5_shockspeed(j,1) = 0;%shockspeedslope(5);
    store_10_shockspeed(j,1) = 0;% shockspeedslope(10);
    store_15_shockspeed(j,1) = 0; %shockspeedslope(15);

    %store_shockspeeds(:,j) = shockspeedslope;
    store_shockspeeds{j,1} = shockspeedslope;
    store_shockspeeds_exp{j,1} = shockspeedslope_exp;

    % store max force
    store_max_force(j,1) = max(F_shocktrack_mean);

    % store height, length
    store_heights{j,1} = 100*height_out; % cm
    store_channels{j,1} = channel/100; % m


    % store radial sensor height vs time
    store_height_radials{j,1} = store_heights_bin;

    % store sensor peaks to plot as one
    store_sensor_force_peak_time(j,1) = sensor_force_peak_time;
    store_sensor_force_peak(j,1) = sensor_force_peak;

    %


end


% plot the expected wave speed based on initial overlap
% the proportionality law from Gomez 2013 and others is c~delta0^(1/4)
% c_propto_delta = init_overlap_mean_vs_depth.^(0.25);
% mass = 4/3 * pi *r_p(j).^3 * rho_p(j);
% c_calc = (sqrt(3)/2*(6*r_p(j).^2)*(E_p(j).*(2*r_p(j)).^(0.5) / (3*(1-0.2^2)*mass)) * init_overlap_mean_vs_depth.^(0.5)).^(0.5);
% set(0, 'CurrentFigure', h_shockdetails_speed);
% hold on
% plot(c_calc(2:(length(depthplot)-2)),depthplot(2:(length(depthplot)-2))/(100*r_p(j)),'b--','LineWidth',2,'DisplayName', 'c_0')

%plot(145*c_propto_delta(2:(length(depthplot)-2)),depthplot(2:(length(depthplot)-2))/(100*r_p(j)),'b--','LineWidth',2,'DisplayName', 'c ~ \delta_0^{1/4}')


% set(0, 'CurrentFigure', h_shockdetails_speed_exp);
% hold on
% plot(105*c_propto_delta,depthplot,'b--','LineWidth',2,'DisplayName', 'c ~ \delta_0^{1/4}')

% plot the dilation bars



%keyboard
% set(0, 'CurrentFigure', h_heightdilation)
% yl = get(gca,'YLim')
% axes('NextPlot','add')
% y_bar_plot = 0.8*yl; %[-0.2,0.4];
% subplot(1,1,1)
% %dil_start = 101;  % 1m
% %dil_end = 301; % 1m, 50cm
% plot(ones(1,2)*channel_old(dil_start)/(r_p_old*100),y_bar_plot,'r--',ones(1,2)*channel_old(dil_end)/(r_p_old*100),...
%     y_bar_plot,'r--','LineWidth',2,'HandleVisibility','off')



% plot arrival and peak times as single entries
% set(0, 'CurrentFigure', h_sample_sensor_force);
% axes('NextPlot','add');
% subplot(1,1,1)
% semilogy(store_sensor_force_peak_time,abs(store_sensor_force_peak),'ko','LineWidth',2','DisplayName','Identified Peak')
% %semilogy(timeplot(store_arrival_ind), store_arrival_val,'ro','LineWidth',2,'DisplayName','First Arrival','MarkerSize',10)
% %semilogy(store_peak_time,store_peak_val,'rx','DisplayName','Peak Force','LineWidth',2,'MarkerSize',10)
% legend('show','location','best')
% set(gca,'yscale','log')



% % plot arrival and peak times as single entries
% set(0, 'CurrentFigure', h_sample_sensor_force);
% axes('NextPlot','add');
% subplot(1,1,1)
% hold on
% % semilogy(timeplot(store_arrival_ind), store_arrival_val,'ro','LineWidth',2,'DisplayName','First Arrival','MarkerSize',10)
% % semilogy(store_peak_time,store_peak_val,'rx','DisplayName','Peak Force','LineWidth',2,'MarkerSize',10)
% legend('show','location','best')
% set(gca,'yscale','log')

% set(0, 'CurrentFigure', h_sample_sensor_ystress);
% axes('NextPlot','add');
% subplot(1,1,1)
% hold on
% %semilogy(timeplot(store_arrival_ind), store_arrival_val,'ro','LineWidth',2,'DisplayName','First Arrival','MarkerSize',10)
% semilogy(store_peakstress_time,store_peakstress_val,'ro','DisplayName','Peak Stress','LineWidth',2,'MarkerSize',10)
% legend('show','location','best')
% legend('off')
% set(gca,'yscale','log')

% find the variance in the steady height dilations, computed through height
% change

% znumgrid = 200% length(height_pdiff_out);
% pdiff_height_mean = zeros(1,znumgrid);
% pdiff_height_variance = zeros(1,znumgrid);
% for jj = 1:znumgrid
%     temp_pdiff_height = [store_ssdilation_heights{1,1}(jj), ...
%             store_ssdilation_heights{2,1}(jj), ...
%             store_ssdilation_heights{3,1}(jj)];
%
%         pdiff_height_variance(1,jj) = std(temp_pdiff_height);
%         pdiff_height_mean(1,jj) = mean(temp_pdiff_height);
% end
% h_pdiff_height_var = figure;
% errorbar(channel(1:znumgrid)/100,pdiff_height_mean,pdiff_height_variance,'bo','DisplayName','Mean Height')
% grid on
% hold on
% set(gca,'FontSize',14)
% ylabel('Mean Percent Change in z')
% xlabel('R [m]')






% call function to produce volume fraction error plots
% [depths_dil,depth_std_cell, channels,radial_std_cell, radial_std_cell_xtop, ...
%     h_errbydepth,h_errbyradial,h_errbyradial_xtop] = ...
%     ErrFun_vf(store_ssdilation, store_channel, store_depthplot,comptype);
%
% saveas(h_errbydepth, [plotpath 'dilationerrorvsdepth.png'])
% saveas(h_errbyradial, [plotpath 'dilationerrorvsradial.png'])
% saveas(h_errbyradial_xtop, [plotpath 'dilaitonerrorvsradial_xtop.png'])


% call function to produce wave runout plots
% [depths_wave,runout_std,speed_std, h_runouterrbydepth,h_speederrbydepth] = ...
%     ErrFun_wave(store_shockrunout,store_shockspeed, store_depthplot,comptype);
%
% % saveas(h_runouterrbydepth, [plotpath 'waverunouterrorvsdepth.png'])
% saveas(h_speederrbydepth, [plotpath 'wavespeederrorvsradial.png'])


% get the heights and wave speeds vs depth only - turned off at the moment
% if 0
%     [channels,radial_avg,radial_std,mean_err_vf,mean_std_vf] = ErrFun_height(store_heights, store_channels,r_start_ind,r_end_buff);
%
%     [depths,depth_avg,depth_std,mean_err_ws,mean_std_ws] = ErrFun_wavespeed(store_shockspeed, store_depths);
% end
% store_shockspeed{j,1}

% save([harvestpath timestring '_' comptype '_errorvar.mat'], ...
%     'depths_dil','depth_std_cell', 'channels','radial_std_cell','radial_std_cell_xtop', ...
%     'depths_wave','runout_std','speed_std')




% saveas(h_errbydepth, [plotpath 'errorvsdepth.png'])
% saveas(h_errbyradial, [plotpath 'errorvsradial.png'])
% saveas(h_errbyradial_xtop, [plotpath 'errorvsradial_xtop.png'])
%
% save([harvestpath timestring '_' comptype '_errorvar.mat'], 'depths','depth_std_cell', 'channels','radial_std_cell','radial_std_cell_xtop')



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

% normalize by mean lunar grain size
% within the range of 40 to 800 microns as per McKay p 306
r_lunar = 250E-6;

% set(0, 'CurrentFigure', h_dil_vsdepth)
% plot(100*deltarho,(z/100)/r_lunar,'k-','LineWidth',2,'DisplayName','Expected')
% ylim([0,40])

set(0, 'CurrentFigure', h_bulkdilationbinvsdepth)
plot(100*deltarho,(z/100)/r_lunar,'k-','LineWidth',2,'DisplayName','Expected')
ylim([0,170])
current_xl = xlim;
new_xl = [current_xl(1), 5];
xlim(new_xl);
h_thisfig = figure(h_bulkdilationbinvsdepth);
thisfig_pos = h_thisfig.Position;
new_pos = [thisfig_pos(1:2),0.5*thisfig_pos(3), thisfig_pos(4)];
set(gcf,'position',new_pos);
legend('show','location','best')
figure(h_bulkdilationbinvsdepth)
%xlim([-10,4])
%keyboard
%
% % new figure that plots our densitites vs Bandfield densities
% rho_p = 2.5; % enter particle bulk density that was input into sim g/cc
% h_density_vs_depth = figure;
% plot(d_background/1000,z,'k-','LineWidth',2,'DisplayName','Background')
% hold on
% grid on
% plot(d_coldspot/1000,z,'k--','LineWidth',2,'DisplayName','LCS')
% plot(rho_p*vf_rowavg_init,depthplot,'bo','LineWidth',2,'DisplayName','\rho_0')
% plot(rho_p*vf_rowavg_final,depthplot,'b+','LineWidth',2,'DisplayName','\rho_f')
% xlabel('Bulk density [g/cm^3]')
% ylabel('Depth [cm]')
% set(gca,'Ydir','reverse')
% %title(['Depth vs Steady state dilation - ' comptype])
% set(gca,'fontsize',14)
% legend('show','location','best')
%
% set(0, 'CurrentFigure', h_density_vs_depth)
% plot(d_coldspot/1000,(z/100)/r_lunar,'k--','LineWidth',2,'DisplayName','LCS')
% plot(d_background/1000,(z/100)/r_lunar,'k-','LineWidth',2,'DisplayName','Background')
% set(gca,'Ydir','reverse')
% %title(['Depth vs Steady state dilation - ' comptype])
% set(gca,'fontsize',14)
% legend('show','location','best')


% find the mean and std dev dilations
dil_mean = 100*mean(deltarho)*ones(1,length(channel));
dil_std = 100*std(deltarho)*ones(1,length(channel));

% set(0, 'CurrentFigure', h_bulkdilationvsradial)
% plot(channel./100,dil_mean,'k-','LineWidth',2,'DisplayName','Mean CS Dilation')
% plot(channel./100,dil_mean + dil_std,'k--','LineWidth',2,'DisplayName','Std.Dev. CS Dilation')
% plot(channel./100,dil_mean - dil_std,'k--','LineWidth',2,'DisplayName','Std.Dev. CS Dilation')

yplot_lnl = store_avg_shockspeed/store_avg_shockspeed(vel==1)
% linear nonlinear plot
ref2 = vel.^(1/2);
ref4 = vel.^(1/4);
ref5 = vel.^(1/5); % power law - relative
ref6 = vel.^(1/6);
% h_sonicv_rel = figure;
% errorbar(evalvals,store_avg_shockspeed/store_avg_shockspeed(evalvals==1),store_std_shockspeed/store_avg_shockspeed(evalvals==1),'bo--','LineWidth',2,'DisplayName','Peak Force')
% grid on
% hold on
% errorbar(evalvals,store_avg_shockspeed_arr/store_avg_shockspeed_arr(evalvals==1),store_std_shockspeed_arr/store_avg_shockspeed_arr(evalvals==1),'co--','LineWidth',2,'DisplayName','First Arrival')
% plot(evalvals((Lf-5):Lf),ref5((Lf-5):Lf),'k--','LineWidth',2,'DisplayName','c_w ~ v_{p}^{1/5}')
% plot(evalvals((Lf-5):Lf),ref6((Lf-5):Lf),'k-.','LineWidth',2,'DisplayName','c_w ~ v_{p}^{1/6}')
% xlabel('Impact Speed [m/s]')
% ylabel('Sound Speed [m/s]')
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% legend('show','location','best')


% h_sonicv_raw_exp = figure;
% errorbar(vel,store_avg_shockspeed_exp,store_std_shockspeed_exp,'bo','LineWidth',2,'DisplayName','Peak Stress (YY) Speed')
% grid on
% hold on
% %errorbar(evalvals,store_avg_shockspeed_arr,store_std_shockspeed_arr,'co--','LineWidth',2,'DisplayName','First Arrival')
% %errorbar(evalvals,store_avg_shockspeed_max,store_std_shockspeed_max,'ro--','LineWidth',2,'DisplayName','Max Force')
% %plot(evalvals,ref,'k--','LineWidth',2,'DisplayName','c_w ~ v_{im}^{1/5}')
% if velyesno
% plot(vel((Lf-2):Lf),8*ref5((Lf-2):Lf),'k-.','LineWidth',2,'DisplayName','c_w ~ v_{p}^{1/5}')
% plot(vel((Lf-2):Lf),4*ref4((Lf-2):Lf),'k:','LineWidth',2,'DisplayName','c_w ~ v_{p}^{1/4}')
% end
% %plot(evalvals((Lf-3):Lf),6*ref2((Lf-3):Lf),'k-.','LineWidth',2,'DisplayName','c_w ~ v_{p}^{1/2}')
% xlabel('v_p [m/s]')
% ylabel('c_w [m/s]')
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% legend('show','location','best')


% plot the trend at each depth

if velyesno
    depth_print_range = 3:(ssL-1);
    speedplotcolor = turbo(length(depth_print_range));
    speedcolorcount = 1;
    h_speed_trend_each_depth_exp = figure;
    %skip the first depth
    for jj = depth_print_range

        speeds_to_plot = zeros(1,Lf);
        for kk = 1:Lf
            speeds_to_plot(1,kk) = store_shockspeeds_exp{kk,1}(jj,1);
        end

        axes('NextPlot','add');
        subplot(1,1,1)
    %     plot(evalvals,store_shockspeeds(jj,:),'--o','LineWidth',2,'DisplayName',['Depth = ' num2str(round(depthplot(depth_print_range(speedcolorcount)))) ' cm'],'Color',speedplotcolor(speedcolorcount,:))
        plot(vel,speeds_to_plot,'--o','LineWidth',2,'DisplayName',['Depth = ' num2str(round(depthplot(depth_print_range(speedcolorcount)))) ' cm'],'Color',speedplotcolor(speedcolorcount,:))
        hold on
        grid on
        xlabel('v_p [m/s]')
        ylabel('c_w [m/s]')
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')
        title('exp')

        speedcolorcount = speedcolorcount +1;

    end

    plot(vel((Lf-2):Lf),8*ref5((Lf-2):Lf),'k-.','LineWidth',2,'DisplayName','c_w ~ v_{p}^{1/5}')
    plot(vel((Lf-2):Lf),4*ref4((Lf-2):Lf),'k:','LineWidth',2,'DisplayName','c_w ~ v_{p}^{1/4}')
    %plot(evalvals((Lf-3):Lf),2.5*ref2((Lf-3):Lf),'k-.','LineWidth',2,'DisplayName','c_w ~ v_{p}^{1/2}')
    legend('show','location','best')
    legend('off')
end






% h_sonicv_raw = figure;
% errorbar(vel,store_avg_shockspeed,store_std_shockspeed,'bo','LineWidth',2,'DisplayName','Peak Stress (YY) Speed')
% grid on
% hold on
% %errorbar(evalvals,store_avg_shockspeed_arr,store_std_shockspeed_arr,'co--','LineWidth',2,'DisplayName','First Arrival')
% %errorbar(evalvals,store_avg_shockspeed_max,store_std_shockspeed_max,'ro--','LineWidth',2,'DisplayName','Max Force')
% %plot(evalvals,ref,'k--','LineWidth',2,'DisplayName','c_w ~ v_{im}^{1/5}')
% if velyesno
% plot(vel((Lf-2):Lf),8*ref5((Lf-2):Lf),'k-.','LineWidth',2,'DisplayName','c_w ~ v_{p}^{1/5}')
% plot(vel((Lf-2):Lf),4*ref4((Lf-2):Lf),'k:','LineWidth',2,'DisplayName','c_w ~ v_{p}^{1/4}')
% end
% %plot(evalvals((Lf-3):Lf),6*ref2((Lf-3):Lf),'k-.','LineWidth',2,'DisplayName','c_w ~ v_{p}^{1/2}')
% xlabel('v_p [m/s]')
% ylabel('c_w [m/s]')
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% legend('show','location','best')


% plot the trend at each depth

if velyesno
    depth_print_range = 3:(ssL-1);
    speedplotcolor = turbo(length(depth_print_range));
    speedcolorcount = 1;
    h_speed_trend_each_depth = figure;
    %skip the first depth
    for jj = depth_print_range

        speeds_to_plot = zeros(1,Lf);
        for kk = 1:Lf
            speeds_to_plot(1,kk) = store_shockspeeds{kk,1}(jj,1);
        end

        axes('NextPlot','add');
        subplot(1,1,1)
    %     plot(evalvals,store_shockspeeds(jj,:),'--o','LineWidth',2,'DisplayName',['Depth = ' num2str(round(depthplot(depth_print_range(speedcolorcount)))) ' cm'],'Color',speedplotcolor(speedcolorcount,:))
        plot(vel,speeds_to_plot,'--o','LineWidth',2,'DisplayName',['Depth = ' num2str(round(depthplot(depth_print_range(speedcolorcount)))) ' cm'],'Color',speedplotcolor(speedcolorcount,:))
        hold on
        grid on
        xlabel('v_p [m/s]')
        ylabel('c_w [m/s]')
        set(gca, 'XScale', 'log')
        set(gca, 'YScale', 'log')

        speedcolorcount = speedcolorcount +1;

    end

    plot(vel((Lf-2):Lf),8*ref5((Lf-2):Lf),'k-.','LineWidth',2,'DisplayName','c_w ~ v_{p}^{1/5}')
    plot(vel((Lf-2):Lf),5*ref4((Lf-2):Lf),'k:','LineWidth',2,'DisplayName','c_w ~ v_{p}^{1/4}')
    %plot(evalvals((Lf-3):Lf),2.5*ref2((Lf-3):Lf),'k-.','LineWidth',2,'DisplayName','c_w ~ v_{p}^{1/2}')
    legend('show','location','best')
    legend('off')
end






% keyboard
% % figure for abstract
% h_densityprofile_and_change = figure;
% subplot(1,2,1)
% plot(d_coldspot./1000,z,'k--','LineWidth',2,'DisplayName','Cold Spot')
% hold on
% grid on
% plot(d_background./1000,z,'k-','LineWidth',2,'DisplayName','Background Regolith')
% ylabel('Depth [cm]')
% xlabel('Density [g/cm^3]')
% title('Figure 1a) Density Profile')
% set(gca,'Ydir','reverse')
% legend('show','Location','best')
% set(gca,'fontsize',14)
% subplot(1,2,2)
% plot(100*deltarho,z,'k-','LineWidth',2,'DisplayName','Expected')
% hold on
% grid on
% ylabel('Depth [cm]')
% xlabel('% Change in Bulk Density')
% title('Figure 1b) Expected Dilation')
% set(gca,'Ydir','reverse')
% set(gca,'fontsize',14)
% saveas(h_densityprofile_and_change, [plotpath 'densityprofile_and_change.png'])
%
%
% keyboard
    % and use this at the end to mess with the size
    % initial conditions - vs radial
    %set(h_init_col, 'Units', 'pixels', 'Position', [100, 300, 1400, 400]);
%     set(0, 'CurrentFigure', h_init_forcevsradial)
%     legend('show','Location','best')
%     saveas(h_init_forcevsradial, [plotpath 'IC-forceVSradial.png'])

%     set(0, 'CurrentFigure', h_init_velocityvsradial)
%     legend('show','Location','best')
%     saveas(h_init_velocityvsradial, [plotpath 'IC-velocityVSradial.png'])
%
%     set(0, 'CurrentFigure', h_init_VFvsradial)
%     legend('show','Location','best')
%     saveas(h_init_VFvsradial, [plotpath 'IC-velocityVSradial.png'])

    % initial conditions - vs depth
    %set(h_init_forcevsdepth, 'Units', 'pixels', 'Position', [100, 300, 1400, 400]);
%     set(0, 'CurrentFigure', h_init_forcevsdepth)
%     legend('show','Location','best')
%     saveas(h_init_forcevsdepth, [plotpath 'IC-forceVSdepth.png'])

%     set(0, 'CurrentFigure', h_init_velocityvsdepth)
%     legend('show','Location','best')
%     saveas(h_init_velocityvsdepth, [plotpath 'IC-velocityVSdepth.png'])

    set(0, 'CurrentFigure', h_init_VFvsdepth)
    legend('show','Location','best')
    saveas(h_init_VFvsdepth, [plotpath 'IC-vfVSdepth.png'])

    % dilations
    %set(h_dil_full, 'Units', 'pixels', 'Position', [500, 200, 1000, 500]);
%     set(0, 'CurrentFigure', h_dil_vsdepth)
%     legend('show','Location','best')
%     saveas(h_dil_vsdepth, [plotpath 'dilationsVSdepth.png'])
%
%     set(0, 'CurrentFigure', h_dil_vsradial)
%     legend('show','Location','best')
%     saveas(h_dil_vsradial, [plotpath 'dilationsVSradial.png'])
%
%     set(h_dil_partial, 'Units', 'pixels', 'Position', [500, 200, 1000, 500]);
%     set(0, 'CurrentFigure', h_dil_partial)
%     legend('show','Location','best')
%     saveas(h_dil_partial, [plotpath 'partialdilations.png'])

    % maximums
    %set(h_max_averages,'Units', 'pixels', 'Position', [500, 200, 1200, 600]);
%     set(0, 'CurrentFigure', h_max_force_logref)
%     legend('show','Location','best')
%     saveas(h_max_force_logref, [plotpath 'maxforce_logref.png'])

    set(0, 'CurrentFigure', h_max_force)
    legend('show','Location','best')
    saveas(h_max_force, [plotpath 'maxforce.png'])

%     set(0, 'CurrentFigure', h_max_velocity)
%     legend('show','Location','best')
%     saveas(h_max_velocity, [plotpath 'maxvelocity.png'])


    % wave details
    %set(h_shockdetails_F,'Units', 'pixels', 'Position', [500, 200, 1250, 450]);
%     set(0, 'CurrentFigure', h_shockdetails_force)
%     legend('show','Location','best')
%     saveas(h_shockdetails_force, [plotpath 'shockforce.png'])

%     set(0, 'CurrentFigure', h_shockdetails_runout)
%     legend('show','Location','best')
%     saveas(h_shockdetails_runout, [plotpath 'shockrunout.png'])

%     set(0, 'CurrentFigure', h_shockdetails_speed)
%     legend('show','Location','best')
%     saveas(h_shockdetails_speed, [plotpath 'shockspeed.png'])
%
%     set(0, 'CurrentFigure', h_shockdetails_speed_exp)
%     legend('show','Location','best')
%     saveas(h_shockdetails_speed_exp, [plotpath 'shockspeed_exp.png'])
%
%     set(0, 'CurrentFigure', h_shockdetails_speed_arrival)
%     legend('show','Location','best')
%     saveas(h_shockdetails_speed_arrival, [plotpath 'shockspeed_arrival.png'])
%
%     set(0, 'CurrentFigure', h_shockdetails_speed_maxforce)
%     legend('show','Location','best')
%     saveas(h_shockdetails_speed_maxforce, [plotpath 'shockspeed_maxforce.png'])

    % dilation bins
%     set(h_dilationbins,'Units', 'pixels', 'Position', [500, 200, 1250, 450]);
%     set(0, 'CurrentFigure', h_dilationbins)
%     legend('show','Location','best')
%     saveas(h_dilationbins, [plotpath 'dilationbins.png'])
%keyboard

    set(0, 'CurrentFigure', h_bulkdilationbinvsdepth)
    legend('show','Location','best')
    saveas(h_bulkdilationbinvsdepth, [plotpath 'bulkdilation_vsdepth.png'])

%     set(0, 'CurrentFigure', h_bulkdilationvsradial)
%     legend('show','Location','best')
%     saveas(h_bulkdilationvsradial, [plotpath 'bulkdilation_vsradial.png'])

%% write to file for saving
    datetime.setDefaultFormats('default','yyyy-MM-dd-HH-mm-ss')

    currenttime = datetime;
    formatOut = 'mm_dd_yyyy';
    timestring = datestr(currenttime,'yyyy-mm-dd-HH-MM-ss');

    % directory of write files, needs to be manually changed
    filedir_write = '/Users/hartzell_lab/DOCUMENTS/COLDSPOTS/MATLAB ANALYSIS/ShockInducedDilation/Prod_v/write_files/';

    file_write_name = strcat(comptype,'_',convertCharsToStrings(timestring),'.txt');
    file_write = strcat(filedir_write,file_write_name);
    fileID = fopen(file_write,'wt');

    fprintf(fileID,'#### channel details #### \n');
    printer_fun(fileID,depthplot,'depthplot');
    printer_fun(fileID,channel,'channel');

    fprintf(fileID,'\n#### particle details #### \n');
    printer_fun(fileID,E_p,'E_p');
    printer_fun(fileID,r_p,'r_p');
    printer_fun(fileID,rho_p,'rho_p');

    fprintf(fileID,'\n#### wave speed details #### \n');
    printer_fun(fileID,vel,'vels');
    printer_fun(fileID,store_avg_shockspeed,'avg_shockspeed');
    printer_fun(fileID,store_std_shockspeed,'std_shockspeed');
    printer_fun(fileID,store_avg_shockspeed_exp,'cprop_avg_');
    printer_fun(fileID,store_std_shockspeed_exp,'cprop_std_');

    printer_fun(fileID,store_solitary_overlap_mean,'solitary_overlap_mean');
    printer_fun(fileID,store_solitary_overlap_std,'solitary_overlap_std');


    fprintf(fileID,'\n#### dilation details #### \n');
    printer_fun(fileID,store_ssdilation_height_mean,'store_ssdilation_height_mean');
    printer_fun(fileID,store_ssdilation_height_std,'store_ssdilation_height_std');
    printer_fun(fileID,store_init_overlap_mean,'store_init_overlap_mean');
    printer_fun(fileID,store_init_overlap_std,'store_init_overlap_std');
    printer_fun(fileID,cell2mat(store_initial_avg_phi),'store_initial_avg_phi');
    printer_fun(fileID,cell2mat(store_initial_std_phi),'store_initial_std_phi');
    printer_fun(fileID,store_loft_depth,'store_loft_depth');
    printer_fun(fileID,cell2mat(store_initial_std_phi),'store_initial_std_phi');

    fprintf(fileID,'\n#### force details #### \n');
    % all solitary wave forces
    printer_fun(fileID,store_Fm_avg,'store_Fm_avg');
    printer_fun(fileID,store_Fm_std,'store_Fm_std');
    printer_fun(fileID,store_F_SW_avg,'store_F_SW_avg');
    printer_fun(fileID,store_F_SW_std,'store_F_SW_std');
    printer_fun(fileID,store_Fm_scaled_avg,'store_Fm_scaled_avg');
    printer_fun(fileID,store_Fm_scaled_std,'store_Fm_scaled_std');
    printer_fun(fileID,store_Fm_pal,'store_Fm_pal');
    printer_fun(fileID,store_F_SW_Fc_avg,'store_F_SW_Fc_avg');
    printer_fun(fileID,store_F_SW_Fc_std,'store_F_SW_Fc_std');
    printer_fun(fileID,store_F_SW_Fc_analytic_avg,'store_F_SW_Fc_analytic_avg');
    printer_fun(fileID,store_F_SW_Fc_analytic_std,'store_F_SW_Fc_analytic_std');

    fprintf(fileID,'\n#### particle velocity #### \n');
    printer_fun(fileID,store_solitary_pvelocity_mean,'store_solitary_pvelocity_mean');
    printer_fun(fileID,store_solitary_pvelocity_std,'store_solitary_pvelocity_std');

    % blast loading
    fprintf(fileID,'\n#### blast loading #### \n');
    printer_fun(fileID,store_BL_m,'store_BL_m');
    printer_fun(fileID,store_BL_nump,'store_BL_nump');
    printer_fun(fileID,store_BL_time,'store_BL_time');
    printer_fun(fileID,store_BL_decayspeed,'store_BL_decayspeed');
    printer_fun(fileID,store_BL_particleV,'store_BL_particleV');

%     % decay
%     fprintf(fileID,'\n#### decay #### \n');
%     %printer_fun(fileID,store_normalized_channel,'store_normalized_channel');
%     printer_fun(fileID,store_dm_normalized,'store_dm_normalized');
%     printer_fun(fileID,store_dm_time,'store_dm_time');
%     store_normalized_channel = cell(Lf,1);
%     store_dm_normalized = cell(Lf,1);
%     store_dm_time{j,1} = cell(Lf,1);

    fclose(fileID);



function [drhoz,drho1,drho2] = densityVdepth(rhod,rho0,H1,H2,z)
% See Bandfield et al, 2014 - Coldspots (Icarus)

drho1 = rhod - (rhod - rho0) * exp(-z/H1);
drho2 = rhod - (rhod - rho0) * exp(-z/H2);


drhoz = (drho2)/(drho1) - 1 ;

end
