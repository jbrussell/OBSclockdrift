% Use daily cross-correlations functions to measure clock drift. One of the two
% stations must have correct time, otherwise this will not work. Daily
% drift is determined by cross-correlating daily CCFs with the reference
% stack.
%
% See Hable et al. 2018 (GJI) [doi:10.1093/gji/ggy236]
%
% jbrussell 8/6/2020
clear
setup_parameters;
IsFigure = 0;
IsFigure_GAUS = 0; % Plot frequency domain filtered and unfiltered

%======================= PARAMETERS =======================%
comp = 'PP'; %'ZZ'; %'RR'; %'TT';
coperiod = [3 8]; % Periods to filter between

stas1 = {'WC01','EC01'}; % station without drift correction
stas2 = {'CC04','CC04'}; % reference station (with correct timing)

% Parameters for making reference stack
Ndaystack = 30; % Number of days after start day to consider in reference stack

% Parameters for estimating drift from
Nsmooth = 7+8; % days to smooth over (MUST BE ODD)
SNRthresh = 2; % Only consider traces with SNR >= SNRthresh
CCthresh = 0.2; % Lowest coherence to consider

%%% --- Parameters to build up gaussian filters --- %%% 
% (effects the width of the filter in the frequency domain)
costap_wid = 0.25; % 0 => box filter; 1 => Hann window

isplotwin = 1;
% Window Velocities
max_grv = 5.5; %8.0;
min_grv = 1.0; %2.2; %1.6; % FOR WINDOWING!

isploth20 = 0;
h20_grv = 1.5;
%==========================================================%
xlims = [-500 500];

dt = parameters.dt;
stalist = parameters.stalist;
nsta = length(stalist);
winlength = parameters.winlength;
figpath = parameters.figpath;
fig_winlength_path = figpath;

% %------------ PATH INFORMATION -------------%
ccf_path = parameters.ccfpath;
ccf_winlength_path = [ccf_path,'window',num2str(winlength),'hr/'];
ccf_daystack_path = [ccf_winlength_path,'dayStack/'];

% ccf_stack_path = ccf_fullstack_path;
ccf_daystack_path = ccf_daystack_path;

figpath = [fig_winlength_path,'clock_drift/',num2str(coperiod(1)),'_',num2str(coperiod(2)),'s/'];
% create figure directory
if ~exist(fig_winlength_path)
    mkdir(fig_winlength_path)
end
if ~exist(figpath)
    mkdir(figpath)
end

%------------ LOAD DATA AND PLOT IN TIME DOMAIN -------------%
for ista1= 1:length(stas1) %1%1:nsta % loop over all stations
    sta1 = stas1{ista1};
    sta2 = stas2{ista1};
    npairall = 0;
    ccf_stack = {};
    nstapair = 0;

    % if same station, skip
    if(strcmp(sta1,sta2))
        continue
    end
    nstapair = nstapair + 1;
    npairall = npairall + 1; % number of total station pairs
    
    % Get correct station lats & lons
    Ista1 = find(strcmp(stalist,sta1));
    Ista2 = find(strcmp(stalist,sta2));
    sta1sta2_dist(nstapair) = deg2km(distance(stalat(Ista1),stalon(Ista1),stalat(Ista2),stalon(Ista2)));
    
    %% Load dayfiles and do smartstack
    ccf_day_path = [ccf_daystack_path,'ccf',comp,'/',sta1,'/',sta1,'_',sta2,'_*.mat'];
    fils = dir(ccf_day_path);
    din = [];
    timeflag = [];
    tvec = NaT(length(fils),1);
    daySNR = [];
    for ifil = 1:length(fils)
        fldr = [ccf_daystack_path,'ccf',comp,'/',sta1,'/'];
        temp = load([fldr,'/',fils(ifil).name]);
        tvec(ifil) = temp.starttime;
        ccf_day = temp.coh_sum./temp.coh_num_day;
        [ ccf_filtered ] = tukey_filt( ccf_day,coperiod,dt,costap_wid );

%             % Do butterworth filtering after rearranging
%             [b, a] = butter(4,flip(1./coperiod)*2*dt); % Butterworth Filter
%             ccf_ifft = filtfilt(b,a,real(ifft(ccf_day)));
%             ccf_filtered = fft(ccf_ifft);

        ccf_day = ccf_filtered;
        N = length(ccf_day);
        ccf_day_ifft = real(ifft(2*ccf_day([1:N/2+1]),N)); % inverse FFT to get time domain
        %rearrange and keep values corresponding to lags: -(len-1):+(len-1)
        ccf_day_ifft = [ccf_day_ifft(end-N+2:end)' ; ccf_day_ifft(1:N)'];

        % Build time axis
        N = length(ccf_day_ifft);
        time = ([0:N-1]-floor(N/2))*dt;
        timeall = [time(time<0), time(time>=0)];
        Ikeep = abs(timeall)<=500; % index data to keep
        timeflag(:,ifil) = timeall(Ikeep);
        din(:,ifil) = ccf_day_ifft(Ikeep);

        % Calculate SNR
        r = sta1sta2_dist(nstapair);
        tmin = -r/min_grv; %50;
        tmax = r/min_grv; %150;
        Isignal = timeflag(:,ifil)>=tmin & timeflag(:,ifil)<=tmax;
        signal = din(Isignal,ifil);
        noise = din(~Isignal,ifil);
        daySNR(ifil) = max(abs(signal))/sqrt(mean(noise.^2));
%             daySNR(ifil) = sqrt(mean(signal.^2))/sqrt(mean(noise.^2));
    end
    dayvec = [0; cumsum(seconds(diff(tvec))/60/60/24)]'; %[1:Ndays];
    % Calculate reference stack
    dstack = sum(din(:,1:Ndaystack),2)./max(abs(sum(din(:,1:Ndaystack),2)));
    ccf_stack{npairall} = dstack;

    twin = r/min_grv;
    if twin < 50
        twin = 50;
    end
    Igrvwin = timeflag(:,1)<=twin & timeflag(:,1)>=-twin;
    Igrvwinmat = repmat(Igrvwin,1,size(din,2));
    norm_mat = repmat(max(abs(din.*Igrvwinmat)),size(din,1),1);
    din_norm = din./norm_mat; din_norm(isnan(din_norm))=0;

    f101 = figure(101); clf; hold on;
    set(gcf,'color','w','position',[340   164   572   995]);
    N= length(ccf_stack{npairall});
    time = [-N/2:N/2]*dt; % build lagtime vector for plotting
    indtime = find(abs(time)<=500);
    imagesc(timeflag(:,1)',dayvec,din_norm');

    plot(time(indtime),dstack(indtime)*10-10,'-k','linewidth',1);
    xlim([-250 250]);
    ylim([-20 max(dayvec)]);
    xlabel('Lag time (s)','fontsize',15);
    ylabel('Day of deployment','fontsize',15);
    title([sta1,'-',sta2,' ',num2str(round(r)),' km (',num2str(coperiod(1)),'-',num2str(coperiod(2)),'s)'],'fontweight','bold','fontsize',15);
    colormap(redbluecmap);
    set(gca,'fontsize',15,'linewidth',1.5,'tickDir','out','box','on','layer','top','ydir','reverse')
    caxis([-1 1]);
    save2pdf([figpath,'ccf',comp,'_',sta1,'_',sta2,'.pdf'],f101,100);

    %% Do cross corelation analysis to measure clock drift

    Ndays = size(din,2);
    CCi = [];
    dtau = [];
    SNRind = ones(size(daySNR));
    SNRind(daySNR<SNRthresh) = nan;
    for iday = 1:Ndays
        % moving window smoothing
        ismooth = iday + [-(Nsmooth-1)/2:(Nsmooth-1)/2];
        ismooth = ismooth(ismooth>=1 & ismooth<=Ndays);
        if Nsmooth == 1
            ismooth = iday;
        end
        iSNRmask = repmat(SNRind(ismooth),size(din,1),1);
        dday = nanmean(din(:,ismooth).*iSNRmask,2);

        % Reference trace
        dref = ccf_stack{npairall};

        % Cross-correlation
        [cc,lags] = xcorr(dday,dref,10/dt);
        lags = lags*dt;
        CC = cc/sqrt((dref'*dref)*(dday'*dday));
        [CCi(iday),Ilag] = max(CC);
        dtau(iday) = lags(Ilag);
    end
    ibadCC = CCi<CCthresh;
    dtau(ibadCC) = nan;
    CCi(ibadCC) = nan;

    % Fit line to data
    c = polyfit(dayvec(~ibadCC),dtau(~ibadCC),1);
    y_est = polyval(c,dayvec);
    secperday = c(1);

    %% plot drift
    f102 = figure(102); clf;
    set(gcf,'color','w');

    ax(1) = subplot(2,1,1); box on;
    plot(tvec,dtau,'-ok','linewidth',2); hold on;
    h102 = plot(tvec,y_est,'-r','linewidth',2);
    ylabel('Clock drift (sec)','fontsize',15);
    title([sta1,'-',sta2,' ',num2str(round(r)),' km (',num2str(coperiod(1)),'-',num2str(coperiod(2)),'s)'],'fontweight','bold','fontsize',15);
    set(gca,'fontsize',15,'linewidth',1.5,'layer','top')
    pos = get(gca,'Position');
    if secperday >= 0
        loc = 'southeast';
    else
        loc = 'northeast';
    end
    lh = legend(h102,[num2str(secperday*1000,'%.1f'),' ms/day'],'location',loc,'Box','off');
    lh.Box = 'off';
    % lh.Position = lh
    set(gca,'Position',pos);

    ax(2) = subplot(2,1,2);
    plot(tvec,CCi,'-ok','linewidth',2);
    xlabel('Date','fontsize',15);
    ylabel('Coherence','fontsize',15);
    set(gca,'fontsize',15,'linewidth',1.5,'layer','top')
    xlim(ax(1).XLim);

    save2pdf([figpath,'ccf',comp,'_',sta1,'_',sta2,'_clockdrift.pdf'],f102,100);

    %% Save drift
    pathsave = './clock_drifts/';
    if ~exist(pathsave)
        mkdir(pathsave);
    end
    filename = [pathsave,'clockdrift',comp,'_',sta1,'_',sta2,'_',num2str(coperiod(1)),'_',num2str(coperiod(2)),'s.mat'];
    save(filename,'secperday','c','tvec','dayvec','CCi','dtau');
end % ista1

