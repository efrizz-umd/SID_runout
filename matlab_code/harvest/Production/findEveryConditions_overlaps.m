function [rowavg, colavg] = ...
    findEveryConditions_overlaps(store,ynumgrids,znumgrids,Ltp,timeplot,depthplot)
%FINDINITIALCONDITIONS Summary of this function goes here
%   Detailed explanation goes here

%% initialize storage
%means
rowavg_init = NaN(1,znumgrids);
colavg_init = NaN(1,ynumgrids);

% standard deviations
rowstd_init = NaN(1,znumgrids);
colstd_init = NaN(1,ynumgrids);


% mean mean
rowavg = zeros(1,Ltp);
colavg = zeros(1,Ltp);


% loop over every time
% plot the vs depth evolution of volume fraction
h_vf_vs_depth_vs_time = figure;

keyboard
numdiv = 1;
plotcolor = turbo(round(Ltp/numdiv));
for ii = 1:Ltp

% loop and find by row or by colum
    % the rows
    for j = 1:znumgrids
        tempR = NaN(1,ynumgrids);

        % build the initial conditions across the row
        for k = 1:ynumgrids
            tempR(1,k) = store{j,k}(ii);


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
            tempC(1,k) = store{k,j}(ii);

        end

        % mean and standard deviations
        colavg_init(1,j) = nanmean(tempC);
        colstd_init(1,j) = nanstd(tempC);    

    end


    % zero protection
    %means
    rowavg_init(rowavg_init == 0) = NaN;
    colavg_init(colavg_init == 0) = NaN;
    
    rowavg(1,ii) = mean(rowavg_init);
    colavg(1,ii) = mean(colavg_init);

    % standard deviations
    rowstd_init(rowstd_init == 0) = NaN;
    colstd_init(colstd_init == 0) = NaN;
    
    if ~rem(ii,numdiv) %&& ii > 150
        set(0, 'CurrentFigure', h_vf_vs_depth_vs_time)
        axes('NextPlot','add')
        subplot(1,1,1)
        plot(rowavg_init,depthplot,'o','LineWidth',1.5,'Color',plotcolor(round(ii/numdiv),:))
        ylabel('Depth [cm]')
        xlabel('\delta [m]')
        set(gca,'Ydir','reverse')
        set(gca,'fontsize',14)
        
        ii
        
        keyboard
    end
    
end

figure
semilogy(timeplot,rowavg,'ko')
hold on
grid on
xlabel('Time [s]')
ylabel('\delta [m]')
set(gca,'fontsize',14)




% and plot the final volume fraction

zzz = zeros(1,length(rowavg_init));
figure
errorbar(rowavg_init,depthplot,zzz,zzz,rowstd_init,rowstd_init,'o','LineWidth',2);
hold on
grid on
ylabel('Depth [cm]')
xlabel('\delta [m]')
set(gca,'Ydir','reverse')
%title(['Depth vs Initial volume fraction - ' comptype])
set(gca,'fontsize',14)

keyboard
end

