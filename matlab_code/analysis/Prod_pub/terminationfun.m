% terminationfun.m
% Eric Frizzell, 2024 - efrizz@umd.edu
% https://github.com/efrizz-umd/SID_runout

function [terminationp,termindex,p_f] = terminationfun(x,y,wspacing,depth,icf)

% ************************************************************************
% This function takes a storage array and computes the averages across columns and row
%
% % ----------- output ----------- %
% - terminationp - value at termination, unused
% - termindex - index of termination point
% - p_f - fit parameters
%
% % ----------- intput ----------- %
% - x - radial positions
% - y - quantity to assess for terminatin (generally, force at each radial pos.)
% - wspacing - size of windows (m)
% - depth - dpeth index
% - icf - initial condition force at that depth
%
% ************************************************************************

if 0 %depth == 5
    figure
    semilogy(x,y,'k.','DisplayName',['Force magnitude']) % ['Force magnitude at depth ' num2str(depth) ' cm']
    hold on
    grid on
    set(gca,'fontsize',14)
    xlabel('Channel Position [m]')
    ylabel('Force Magnitude [N]')
    set(gca,'FontSize',14)
    legend('show','location','best')

end


[maxF,~] = max(y);
[mindist,minind] = min(y);
terminationp = x(minind);
%termindex = minind

% find the first time F < initial condition force prior to the minimum
shockset = y(1:minind);
aboveinit = shockset - icf;
[~,crossoverp] = min(abs(aboveinit));
termindex = crossoverp;
% now from this reduced set find the first time that F within 10% of ICF
%tenp = icf*1.1;
%shockrange = (shockset - icf)/icf;
shockrange = shockset./icf;
% find the first instance of the force being within 10 percent of ICF
min2 = find(shockrange < 100.0,1,'first')-1;
%
% find the first instance of the force being 1 percent of the max
min3 = find(shockset < .01*maxF,1,'first')-1;


fw = 50; % size of filter window
y_smooth = smoothdata(y,'movmean',fw);

if min2 < 1
    min2 = 1;
end

if ~isempty(min2)
    %termindex = min2;
    termindex = min2;
end


% %plot it
% semilogy(x(minind),y(minind),'ro','LineWidth',2,'DisplayName','Termination Point')
% semilogy(x(crossoverp),y(crossoverp),'rx','LineWidth',2,'DisplayName','Termination Point with ICF error')
if 0 %depth == 5
    icf_f_plot = 10*ones(1,length(x)) * icf;
    semilogy(x(min2),y(min2),'rs','LineWidth',2,'DisplayName','Termination Point')
    semilogy(x,icf_f_plot,'--r','LineWidth',2,'DisplayName','10x Initial Force')
    set(gca,'fontsize',14)
    %keyboard
end

% if the 10% point is defined, use it for termination point instead
% if ~isempty(min2)
%     termindex = min2;
% end
%keyboard
% legend('show','location','best')
%
%
% keyboard
dy = abs(finite_diff_d1_fun(y,wspacing));
% find the moving mean of the input data
% mean filter
fw = 50; % size of filter window
dy_mean = smoothdata(dy,'movmean',fw);


% split the data into first 2m, last 2m
% x_f = x(1:200);
% dy_f = dy(1:200);
%x_f = x(30:200);
%dy_f = dy_mean(30:200); % use the smoothed data to fit to
%x_end = x(400:end);
%dy_end = dy(400:end);
% mean over the end points
%mean_end = mean(dy_end);
% for straight line plotting
%mean_end_plot = mean_end*ones(1,length(x));

% fit two functions
% myfirstfun = @(param) norm(param(2)*10.^(-param(1)*(x_f -x_f(1))) + param(3) - dy_f,2);
% myfirstfun = @(param) norm(param(2)*param(3).^(-param(1)*(x_f -x_f(1))) - dy_f,2);
myfirstfun = @(param) norm(param(2)*param(3).^(-param(1)*(x_f -x_f(1))) - dy_f,2);

% initial guess as exponential decay

p0_f = [.05,1,3];
% p_f = fminsearch(myfirstfun,p0_f);
p_f = 0;

% dy_fit_f = p_f(2)*10.^(-p_f(1)*(x_f-x_f(1))) + p_f(3);
% dy_fit_f = p_f(2)*p_f(3).^(-p_f(1)*(x_f-x_f(1))); % front portion
%dy_fit_a = p_f(2)*p_f(3).^(-p_f(1)*(x-x_f(1)));

% find the minimum distance between the two lines
% since the data points are all at the same x, just use y
%dist_abs = abs(dy_fit_a - mean_end_plot);
%[mindist,minind] = min(dist_abs);
%terminationp = x(minind);
%termindex = minind;


% plot for checking
% plot absolute value
    %h_termfinding = figure;
%% single figure
% plot(x,dy,'k.','DisplayName','dF/dy Finite Diff. Estimate')
% grid on
% hold on
% plot(x,dy_mean,'m:','DisplayName','Mean Smoothed','LineWidth',3)
% plot(x,dy_fit_a,'b-.','LineWidth',3,'DisplayName','Log Decay Fit')
% plot(x,mean_end_plot,'g--','LineWidth',3,'DisplayName','Mean Floor')
% plot(terminationp,dy_fit_a(minind),'r*','LineWidth',3,'DisplayName','Termination','MarkerSize',10)
% title(['Estimate of dF/dy and shock termination point at depth ' num2str(depth) 'cm'])
% set(gca,'yscale','log')
% set(gca,'fontsize',14)
% legend('show','location','best')
% ylabel('abs(dF/dy)')
% xlabel('Channel Position [cm]')

%% subplot
% subplot(1,2,1);
% plot(x./100,dy,'k.','DisplayName','dF/dy Finite Diff. Estimate')
% grid on
% hold on
% plot(x./100,dy_mean,'m:','DisplayName','Mean Smoothed','LineWidth',3)
% plot(x./100,dy_fit_a,'b-.','LineWidth',3,'DisplayName','Log Decay Fit')
% plot(x./100,mean_end_plot,'g--','LineWidth',3,'DisplayName','Mean Floor')
% plot(terminationp./100,dy_fit_a(minind),'r*','LineWidth',3,'DisplayName','Termination','MarkerSize',10)
% %title(['Estimate of dF/dy and shock termination point at depth ' num2str(depth) 'cm'])
% set(gca,'yscale','log')
% set(gca,'fontsize',14)
% legend('show','location','best')
% ylabel('abs(dF/dy) [N/m]')
% xlabel('Channel Position [m]')

%  keyboard


end
