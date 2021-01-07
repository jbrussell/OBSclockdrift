% Read .mat drifts and write them to a text file
%
% jbrussell 8/6/2020
clear
setup_parameters;
IsFigure = 0;
IsFigure_GAUS = 0; % Plot frequency domain filtered and unfiltered

%======================= PARAMETERS =======================%
comp = 'PP';
coperiod = [3 8]; % Periods to filter between

stas1 = {'WC01','EC01'}; % station without drift correction
stas2 = {'CC04','CC04'}; % reference station (with correct timing)

% Parameters for making reference stack
Ndaystack = 30; % Number of days after start day to consider in stack

% Parameters for estimating drift from
Nsmooth = 7+8; % days to smooth over (MUST BE ODD)
SNRthresh = 2; % Only consider traces with SNR >= SNRthresh
CCthresh = 0.2; % Lowest coherence to consider


fid = fopen('YoungORCA_drifts.txt','w');
fprintf(fid,'STA       DRIFT (MS/DAY)\n');

figpath = [parameters.figpath(1:end-1),'_correctclock/'];
figpath = [figpath,'clock_drift/',num2str(coperiod(1)),'_',num2str(coperiod(2)),'s/'];

%------------ LOAD DATA AND PLOT IN TIME DOMAIN -------------%
f70 = figure(70); clf;
set(gcf,'Position',[183     9   782   696]);
lgd = {};
clrs = brewermap(length(stas1),'set1');
for ista1= 1:length(stas1) %1%1:nsta % loop over all stations
    sta1 = stas1{ista1};
    sta2 = stas2{ista1};

    %% Load drifts
    pathsave = './clock_drifts/';
    filename = [pathsave,'clockdrift',comp,'_',sta1,'_',sta2,'_',num2str(coperiod(1)),'_',num2str(coperiod(2)),'s.mat'];
    load(filename,'secperday','c','tvec','dayvec','CCi','dtau');
    temp = load(filename);
    
    fprintf(fid,'%s      %f\n',sta1,secperday*1000);
    
    ax(1) = subplot(2,1,1);
    h1(ista1) = plot(temp.tvec,temp.dtau,'-o','color',clrs(ista1,:),'markerfacecolor',clrs(ista1,:),'linewidth',1); hold on;
%     plot(temp.tvec,y_est,'-','color',clrs(ista1,:));
    ylabel('Clock drift (sec)','fontsize',15);
    xlabel('Date');
%     title([num2str(coperiod(1)),'-',num2str(coperiod(2)),'s)'],'fontweight','bold','fontsize',15);
    set(gca,'fontsize',15,'linewidth',1.5,'layer','top')
    
    ax(2) = subplot(2,1,2);
    y_est = polyval(c,temp.dayvec);
    h2(ista1) = plot(temp.tvec,temp.dtau - y_est,'-o','color',clrs(ista1,:),'markerfacecolor',clrs(ista1,:),'linewidth',1); hold on;
%     plot(temp.tvec,y_est,'-','color',clrs(ista1,:));
    ylabel('Residuals from linear fit (sec)','fontsize',15);
    xlabel('Date');
%     title([num2str(coperiod(1)),'-',num2str(coperiod(2)),'s)'],'fontweight','bold','fontsize',15);
    set(gca,'fontsize',15,'linewidth',1.5,'layer','top')
    
    lgd1{ista1} = [sta1,': ',num2str(secperday*1000,'%.1f'),' ms/day'];
    lgd2{ista1} = [sta1,': ',num2str(sqrt(nanmean((temp.dtau - y_est).^2))*1000,'%.1f'),' ms'];
end % ista1
fclose(fid);

plot(ax(1),ax(1).XLim,[0 0],'--k','linewidth',1.5);
plot(ax(2),ax(2).XLim,[0 0],'--k','linewidth',1.5);
lh = legend(h1,lgd1,'location','eastoutside','Box','off');
lh.Box = 'off';
legend(h2,lgd2,'location','eastoutside','Box','off');
ax(1).Position(3) = ax(2).Position(3);
% drawnow;

save2pdf([figpath,'alldrifts_summary.pdf'],f70,100);

