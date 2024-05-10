% least_fit_fun_nearestneighbor.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [time_uniform,gen_uniform,derivative,p] = least_fit_fun_nearestneighbor(time_data,pos_data,evalval,depthind)

% ************************************************************************
% This function takes a storage array and computes the averages across columns and row
%
% % ----------- output ----------- %
% - time_uniform - smooth time array for plotting
% - gen_uniform - smotth dervative for plotting
% - derivative - computed slope
% - p - fit parameters
% % ----------- intput ----------- %
% - time_data - time array (s) (independent variable)
% - pos_data - position array (m) (dependent variable)
% - evalval - sometimes used as a flag, currently commented out
% - depthind - used for plotting, indicates sensor depth
%
% ************************************************************************

% non linear least squares fit

% t_1 = time_data(1);
% p_1 = pos_data(1);
%
% % set cooridinates of power function to origin
% time_data = time_data - t_1;
% pos_data = pos_data - p_1;

mean(time_data);


% myfun = @(x,xdata) x(1)*(xdata).^(x(2));
% p0 = [1,(3/6)]; % initial guess
% [p,~,residuals] = lsqcurvefit(myfun,p0,time_data,pos_data);

% % use evalval <= 1 for compressive tests. Checking for shear
% if  1 %evalval <= 1 %0.2% 1 % 0.2 % reset to linear p(2) >= 1 ||p(2) >= 1 ||
%     % reset to original data
% %     time_data = time_data + t_1;
% %     pos_data = pos_data + p_1;

    % find the average distance between each point and the next
    nump = length(time_data);
    dist_arr = zeros(1,nump);

    % end points have one distance value, otherwise take the average
    dist_arr(1,1) = norm( (time_data(2) - time_data(1)), (pos_data(2) - pos_data(1)) );
    dist_arr(1,nump) = norm( (time_data(nump) - time_data(nump-1)), (pos_data(nump) - pos_data(nump-1)) );



    for jj = 2:(nump-1)
        dist_arr(1,jj) = mean( ...
            [(norm( [(time_data(jj) - time_data(jj-1)), (pos_data(jj) - pos_data(jj-1))] )  ), ...
            (norm( [(time_data(jj+1) - time_data(jj)), (pos_data(jj+1) - pos_data(jj) )] ) )] );

    end

    % 0.05 heuristically chosen
    mean_dist_arr = mean(dist_arr(2:(end-1)));
    stdv_dist_arr = std(dist_arr(2:(end-1)));
    mean_minus_dist_arr = dist_arr - mean_dist_arr;

    dist_nottoobig =  abs(mean_minus_dist_arr) < 2 * stdv_dist_arr;

    % exclude endpoints
%     dist_nottoobig = [0, nottoobig_partial, 0]
%     keyboard
%     dist_nottoobig = dist_arr < 0.05;
    include_time = time_data; %time_data(dist_nottoobig);
    include_pos = pos_data; %pos_data(dist_nottoobig)

    myfun = @(x,xdata) x(1) + x(2)*xdata;
    p0 = [1,5]; % initial guess
%     [p,~,residuals] = lsqcurvefit(myfun,p0,time_data,pos_data);
    [p,~,residuals] = lsqcurvefit(myfun,p0,include_time,include_pos);

    % remove residuals
    good_residuals_ind = abs(residuals) <= 2 * std(residuals);
    bad_residuals_ind = abs(residuals) > 2 * std(residuals);

    good_times = include_time(good_residuals_ind);
    good_pos =  include_pos(good_residuals_ind);


    % new fit with only good points
%     [p,~,~] = lsqcurvefit(myfun,p0,time_data(good_residuals_ind),pos_data(good_residuals_ind));

    [p,~,~] = lsqcurvefit(myfun,p0,good_times,good_pos);


    pol = polyfit(good_times,good_pos,1);
    %keyboard
    p = [pol(2), pol(1)];
%     if isempty(p)
%         [p,~,~] = lsqcurvefit(myfun,p0,good_times,good_pos);
%     end
    dt = 0.00001;
    time_uniform = 0:dt:max(time_data);
    gen_uniform = myfun(p,time_uniform);




% else
%     % remove residuals
%     % find those with points within 2 sigma of the fit
%     good_residuals_ind = abs(residuals) <= 2 * std(residuals);
%     bad_residuals_ind = abs(residuals) > 2 * std(residuals);
%
%     % new fit with only good points
%     [p,~,~] = lsqcurvefit(myfun,p0,time_data(good_residuals_ind),pos_data(good_residuals_ind));
%
%     time_uniform = 0:0.001:max(time_data);
%     gen_uniform = myfun(p,time_uniform);
%
% %     time_uniform = time_uniform + t_1;
% %     gen_uniform = gen_uniform +p_1;
%
% end


%get the derivative
%timestep is 0.001 s
derivative = gradient(gen_uniform,dt);

if isempty(derivative)
    keyboard
end

% % % %% position data
% figure
% plot(time_data,pos_data,'bo','LineWidth',2,'DisplayName',['c ~ x^{' num2str(p(2)) '}']);
% hold on
% grid on
% plot(time_uniform,gen_uniform,'k--','LineWidth',2,'DisplayName','NLLS Fit');
% plot(time_data(bad_residuals_ind),pos_data(bad_residuals_ind),'rx','LineWidth',2,'DisplayName','> 2\sigma');
% xlabel('Time [s]')
% ylabel('Position [m]')
% set(gca,'fontsize',14)
% title(['v_{p} = ' num2str(evalval) 'ms, dpethind ' num2str(depthind)])
% legend('show','location','best')

% figure
% plot(pos_data,dist_arr,'ko')
%


%keyboard
%
% %% wave speed
% %figure
% yyaxis right
% plot(time_uniform,derivative,'ko','LineWidth',2,'DisplayName',['v_{p} = ' num2str(evalval) 'ms, dpethind ' num2str(depthind)]);
% hold on
% grid on
% ylim([0,20]);
% %xlabel('Time [s]')
% ylabel('Speed [m/s]')
% %set(gca,'fontsize',14)
% legend('show','location','best')
%
% %
%  keyboard
% end
