    %% Script used for a Sequential Importance Sampling with Resampling (SISR) Particle Filter using gps data
    %% State Model consists of 5 elements, Measurement Model consists of 2 elements
    %% Generates estimates of each 5 state elements saved in variable xPHatMean
    %% Plots observed/measured states, particle trajectories, particle weights, and the estimated path using the first two state elements (longitude and latitude)

clear all; close all;

rng('default');
name                  = 'Project_5'; %  ADDED

%% LOAD DATA

% LEIF TRIAL AND RUN IN METERS
load('Project/python/leifMeters.mat')
trailStart = [lon_m(1),lat_m(1)];
lon_leif_m = lon_m - trailStart(1);
lat_leif_m = lat_m - trailStart(2);
alt_leif_m = alt_m;

load('Project/python/leifRunAllMeters.mat')
lon_run_m = lon_m - trailStart(1);
lat_run_m = lat_m - trailStart(2);
alt_run_m = alt_m;
time_run = time;
clear lon_m
clear lat_m
clear alt_m
clear time

load('Project/AngleInterpolant_Nearest.mat')

%% for testing

t = 650;
lon_run_m = lon_run_m(450:t); % (800:920);%
lat_run_m = lat_run_m(450:t); % (800:920);%
alt_run_m = alt_run_m(450:t); % (800:920);%
time_run = time_run(450:t-1); % (800:919);%

%% end testing
%% SET UP
SampleTimes = [1 1+cumsum(time_run)];
TotalTime = 1+sum(time_run);

nSamples              = TotalTime;
nParticles            = 5000;
nElements             = 5;

% p(v) based 5k paces for US runners
ps = makedist('stable','alpha',0.75,'beta',-0.85,'gam',2.05,'delta',1.55);
pst = truncate(ps,0,12.5);

% noise variances
t0_sigma = 0.1;
l_sigma = 2; % BASED ON WIDTH OF TRAIL ?!?!?!

% t_c = 1/4;
t_sigma = 0.75; % BASED ON INTERPOLATION ?!?!

v_sigma = 0.9; % velocity sigma
m_sigma = 9; % measurement noise
 
alpha = 0.95;

systemParameters = alpha;
measurementParameters = [];

nEffectiveThreshold = 0.10*nParticles; 

fprintf('Particle:                      %5.3f\n',nParticles);
fprintf('Effective Threshold:           %5.3f\n\n',nEffectiveThreshold);

%% OBSERVATIONS
x = zeros(2,length(lon_leif_m));
y = NaN(2,nSamples);

x(1,:) = lon_leif_m;
x(2,:) = lat_leif_m;
y = NaN(2,nSamples);
y(:,SampleTimes) = [lon_run_m ; lat_run_m];

xMax = 5+max(x(1,:)); 
xMin = -5+min(x(1,:));
yMax = 5+max(y); 
yMin = -5+min(y);

%% PLOT OBSERVATIONS
figure;
subplot(2,1,1)
plot(x(1,:),x(2,:),'b','LineWidth',2);
ylabel('$x_n$','Interpreter','latex');
title('True Position');

subplot(2,1,2)
idxs = find(~isnan(y(1,:)));
plot(y(1,idxs),y(2,idxs),'r-*','LineWidth',2);
ylabel('$y_n$','Interpreter','latex');
title('Observed Position');

%% PARTICLES
datetime
t_start = cputime;
fprintf('particle man, particle man...\n')
% xp = zeros(nParticles,nSamples);
% wp = zeros(nParticles,nSamples);
% xPHatMean = zeros(nSamples,1);

% INITIAL STATE
lon0 = y(1,1); %0;
lat0 = y(2,1); %0;
t0 = A(lon0,lat0); % 2.7433 % +randn*t0_sigma; % -pi+2*pi*rand;
v0 = random(pst);
m0 = random(pst);
initialState = [lon0;lat0;t0;v0;m0];

t = 1;
for n=1:nSamples

    if isnan(y(1,n))
        t = t+1;
    else
        t = 1;
    end

    l_noise = mvnrnd([0,0],[sqrt(t)*l_sigma, 0;0 sqrt(t)*l_sigma],nParticles)';
%     for m = 1:nParticles
%         t_noise(m) = wrappedCauchy(t_c);
%     end
    t_noise = sqrt(t_sigma)*randn(1,nParticles); % random(p,1,nParticles); %
    v_noise = sqrt(v_sigma)*randn(1,nParticles);
    processNoise = [l_noise;t_noise;v_noise;zeros(1,nParticles)];
    
    for cp=1:nParticles
        if n==1
            xp(cp,:,n)      = initialState + processNoise(:,cp);
            xp(cp,3,n)      = wrapTo2Pi(xp(cp,3,n));
            xp(cp,4,n)      = max(xp(cp,4,n),0);
            wp(cp,n)     = 1/nParticles; % TIME 1 = UNIFORM DISTRIBUTION FOR WEIGHTS
        else
            xp(cp,:,n)   = SystemModel(xp(cp,:,n-1),systemParameters,A,n-1)+processNoise(:,cp);
            xp(cp,3,n)   = wrapTo2Pi(xp(cp,3,n));
            xp(cp,4,n)   = max(xp(cp,4,n),0);
            mu(:,:)      = xp(cp,1:2,n);
            if isnan(y(1,n))
                wp(cp,n) = 1/nParticles;
            else
                likelihood   = mvnpdf(y(:,n)',mu,[m_sigma,0;0,m_sigma]);    
                wp(cp,n)     = wp(cp,n-1)*likelihood; % ASSUMES IMPORTANCE DENSITY = PRIOR
            end 
        end
    end            
    wp(:,n) = wp(:,n)/sum(wp(:,n));   
    nEffective(n) = 1/sum(wp(:,n).^2);
    if nEffective(n) < nEffectiveThreshold
        iPick = 1;
        c = cumsum(wp(:,n));
        threshold = rand/nParticles;
        for cp=1:nParticles
            while threshold > c(iPick) && iPick < nParticles, iPick = iPick+1; end;
            xpNew(cp,:) = xp(iPick,:,n);
            threshold = threshold + 1/nParticles;
        end
        xp(:,:,n) = xpNew;
        wp(:,n) = ones(nParticles,1)/nParticles;
    end
    xPHatMean(1,n) = sum(wp(:,n).*xp(:,1,n));
    xPHatMean(2,n) = sum(wp(:,n).*xp(:,2,n));
    xPHatMean(3,n) = sum(wp(:,n).*xp(:,3,n));
    xPHatMean(4,n) = sum(wp(:,n).*xp(:,4,n));
    xPHatMean(5,n) = sum(wp(:,n).*xp(:,5,n));
%     xPHatMean = xPHatMean';
end

t_end = round((cputime - t_start)/60,3);
fprintf('doin the things a particle can...%.2f\n',t_end);

%% EROOR - NEED TO ONLY MEASUERED FOR TIME STAMPS WHERE HAVE MEASUREMENT! 
y_xPHatMean = xPHatMean(1:2,SampleTimes);
y_samples = y(:,SampleTimes);
for r = 1:length(SampleTimes) 
    meas_error(r) = mean((y_samples(1,1:r)-y_xPHatMean(1,1:r)).^2 + (y_samples(2,1:r)-y_xPHatMean(2,1:r)).^2);  
end

%% PLOT PARTICLES
figure;
LonP(:,:) = xp(:,1,:);
plot(LonP','g');%,'Marker','s','MarkerFaceColor','black')
hold on; plot(idxs,y(1,idxs),'r-','Marker','.');
plot(xPHatMean(1,:),'k--');
title('Longitude');

figure;
LatP(:,:) = xp(:,2,:);
plot(LatP','g')
hold on; plot(idxs,y(2,idxs),'r-','Marker','.');
plot(xPHatMean(2,:),'k--');
title('Latitude');

% figure;
% tP(:,:) = xp(:,3,:);
% plot(tP','g')
% diff_lon = diff(lon_run_m);
% diff_lat = diff(lat_run_m);
% dist = sqrt(diff_lon.^2+diff_lat.^2);
% for n = 1 : length(diff_lon)
%     a(n) = atan2(diff_lat(n),diff_lon(n));
% end
% a = wrapTo2Pi(a);
% hold on; plot(idxs(1:end-1),a,'r-','Marker','.');
% plot(xPHatMean(3,:),'k--');
% title('Theta');
% 
figure;
vP(:,:) = xp(:,4,:);
plot(vP','g')
hold on;
distances = sqrt(diff(lon_run_m).^2+diff(lat_run_m).^2);
plot(idxs(1:end-1),distances./diff(SampleTimes),'r-','Marker','.');
plot(xPHatMean(4,:),'k--');
title('Velocity');
% 
% figure;
% mP(:,:) = xp(:,5,:);
% plot(mP','g')
% hold on; 
% % plot(x(5,:),'r');
% plot(xPHatMean(5,:),'k--');
% title('Average Velocity');
% 


%% PLOT WEIGHTS
% figure;
% 
% subplot(3,1,1);
% plot(0:nSamples-1,std(wp),'r');
% ylabel('$\sigma_w$','Interpreter','LaTeX');
% xlim([0 nSamples-1]);
% %ylim([0.1^(nSamples) 1]);
% box off;
% 
% subplot(3,1,2);
% semilogy(0:nSamples-1,wp,'b');
% ylabel('$w_n^i$','Interpreter','LaTeX');
% xlim([0 nSamples-1]);
% ylim([0.1^(nSamples) 1]);
% box off;
%         
% subplot(3,1,3);
% plot(0:nSamples-1,wp,'b');
% xlabel('$n$','Interpreter','LaTeX');
% ylabel('$w_n^i$','Interpreter','LaTeX');
% xlim([0 nSamples-1]);
% ylim([0 1]);
% box off;
% 
%% PLOT ESTIMATE VS MEASURED PATH

figure;hold on;
plot(y(1,:),y(2,:),'r-','Marker','.');
plot(xPHatMean(1,:),xPHatMean(2,:),'g-','Marker','.');
dx = 2; dy=2;
div = 10^(length(num2str(length(y)))-3);
often = round(length(y)/div,0);
for m = 0:often:length(y)
    [idx,l] =  min(abs(SampleTimes-m));
    logic =  (SampleTimes(l)-m)<0;
    if logic
        idx = m-min(abs(SampleTimes-m));
    else
        idx = m+min(abs(SampleTimes-m));
    end
    plot(y(1,idx),y(2,idx),'kp','MarkerSize',10,'MarkerFaceColor','black');
    s = ['n= ',num2str(m)];
    text(y(1,idx)+dx,y(2,idx)+dy,s);
end
legend('Measured','Estimated')
% PlotMovie


%% STATE SPACE MODEL

function [xNext] = SystemModel(x,parameters,A,n)
    a = parameters(1);
    b = binornd(1,0.9);
    x1 = x(1) + (x(4)*cos(x(3)));
    x2 = x(2) + (x(4)*sin(x(3)));
    x3 = x(3); %A(x(1),x(2)); %  
    x4 = b*x(5);
    x5 = a*x(5) + (1-a)*x(4);
    xNext = [x1, x2, x3, x4, x5]';

end

function [y] = MeasurementModel(x,parameters,n)
    y = [x(1), x(2)]';
end

function [b] = KernelFunction(x)
    b = (1-x.^2).^2.*(abs(x)<=1);
end