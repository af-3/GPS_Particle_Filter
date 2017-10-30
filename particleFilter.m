% function [error_dist error_t error_v error_y] = project_4;

    %% Used for Monte Carlo Simulation for a Sequential Importance Sampling with Resampling (SISR) Particle Filter using synthetic data
    %% State Model consists of 5 elements, Measurement Model consists of 2 elements
    
    %% Input:
    
    %% No input required

    %% Output:
    
    %% error_dist = Mean Squared Error of Euclidean Distance between estimated state and true state
    %% error_t = Mean Squared Error of heading (theta) between estimated state and true state
    %% error_v = Mean Squared Error of velocity (v) between estimated state and true state
    %% error_y = Mean Squared Error of Euclidean Distance between measured state and true state, used to evaluate the relative improvement of the estimate



    %% SETUP 
    clear all; close all;


    rng('default');
    name                  = 'Project_4'; %  ADDED

    % Angles for Leif - Interpolated
    load('Project/AngleInterpolant_Nearest.mat','A')
    
    nSamples              = 100;
    nParticles            = 1000;
    nElements             = 5;
    
    nRegions              = 2000;
    kernelWidth           = 0.2;

    % p(v) based 5k paces for US runners
    ps = makedist('stable','alpha',0.75,'beta',-0.85,'gam',2.05,'delta',1.55);
    pst = truncate(ps,0,12.5);
    
    t0_sigma = 0.1;
    l0_sigma = 2;
    l0_noise = mvnrnd([0,0],[l0_sigma,0;0,l0_sigma],1,2);
    lon0 = 0 + l0_noise(1); % beginning of leif @ thurman 6.3508e+06
    lat0 = 0 + l0_noise(2); % 5.0449e+06
    t0 = A(lon0,lat0) + randn*t0_sigma; % based on interopolated angle data
    v0 = random(pst); 
    m0 = random(pst);

    initialState = [lon0;lat0;t0;v0;m0];
    l_sigma = 2; % ~ width of trail

    % t_c = 1/4;
    t_sigma = 0.75; % based on turn angles on leif
    % load('Project/PhiDist.mat');

    v_sigma = 0.9; % velocity sigma

    m_sigma = 9; % measurement noise

    alpha = 0.95;
    systemParameters = alpha;
    measurementParameters = [];

    nEffectiveThreshold = 0.10*nParticles; 

%     fprintf('Particle:                      %5.3f\n',nParticles);
%     fprintf('Effective Threshold:           %5.3f\n',nEffectiveThreshold);

    %% OBSERVATIONS
    x = zeros(nElements,nSamples+1);
    y = zeros(2,nSamples);

    measurementNoise = mvnrnd([0,0],[m_sigma,0;0,m_sigma],nSamples)';

    l_noise = mvnrnd([0,0],[l_sigma, 0;0 l_sigma],nSamples)';
    % for n = 1:nSamples
    %     t_noise(n) = wrappedCauchy(t_c);
    % end
    % t_noise = t_sigma*randn(1,nSamples);
    t_noise = sqrt(t_sigma)*randn(1,nSamples);
    v_noise = sqrt(v_sigma)*randn(1,nSamples);

    processNoise = [l_noise;t_noise;v_noise;zeros(1,nSamples)];

    % CREATE OBSERVATIONS
    x(:,1) = initialState + processNoise(:,1); % ADD IN INITIAL NOISE
    for n=1:nSamples
        x(:,n+1) = SystemModel(x(:,n),systemParameters,A,n) + processNoise(:,n);
        y(:,n  ) = MeasurementModel(x(:,n)',measurementParameters,n)  + measurementNoise(:,n);
    end
    x = x(:,1:nSamples);
    x(3,:) = wrapTo2Pi(x(3,:));
    x(4,:)   = max(x(4,:),0);
    % y(:,10:15) = NaN;

    xMax = 5+max(x(1,:)); 
    xMin = -5+min(x(1,:));
    yMax = 5+max(y); 
    yMin = -5+min(y);

    %% PLOT OBSERVATIONS
    
%     figure;
%     subplot(2,1,1)
%     plot(x(1,:),x(2,:),'b','LineWidth',2);
%     ylabel('Latitude (m)','Interpreter','latex','FontSize',16);
%     xlabel('Longitude (m)','Interpreter','latex','FontSize',16);
%     title('True Position ($x_n$)','Interpreter','latex','FontSize',16);
% 
%     subplot(2,1,2)
%     plot(y(1,:),y(2,:),'r','LineWidth',2);
%     ylabel('Latitude (m)','Interpreter','latex','FontSize',16);
%     xlabel('Longitude (m)','Interpreter','latex','FontSize',16);
%     title('Observed Position ($y_n$)','Interpreter','latex','FontSize',16);

    %% PARTICLES

    t_start = cputime;
    fprintf('Particle man, particle man ...\n');
    % xp = zeros(nParticles,nSamples);
    % wp = zeros(nParticles,nSamples);
    % xPHatMean = zeros(nSamples,1);
    for n=1:nSamples

        l_noise = mvnrnd([0,0],[l_sigma, 0;0 l_sigma],nParticles)';
        t_noise = sqrt(t_sigma)*randn(1,nParticles);
        v_noise = sqrt(v_sigma)*randn(1,nParticles);
        processNoise = [l_noise;t_noise;v_noise;zeros(1,nParticles)];

        for cp=1:nParticles
            if n==1
                xp(cp,1:2,n) = [0,0] + mvnrnd([0,0],[l_sigma,0;0,l_sigma]);
                xp(cp,3,n)   = A(0,0) + randn*t0_sigma; %-pi+2*pi*rand; %initialState + [randn*p_sigma; randn*v_sigma];
                xp(cp,4,n)   = random(pst);
                xp(cp,5,n)   = random(pst);
                wp(cp,n)     = 1/nParticles; % TIME 1 = UNIFORM DISTRIBUTION FOR WEIGHTS
            else
                xp(cp,:,n)   = SystemModel(xp(cp,:,n-1),systemParameters,A,n-1)+processNoise(:,cp);
                xp(cp,3,n)   = wrapTo2Pi(xp(cp,3,n));
                xp(cp,4,n)   = max(xp(cp,4,n),0);
                mu(:,:)      = xp(cp,1:2,n);
                if isnan(y(1,n))
                    wp(cp,n) = 1/nParticles;
                else
                    likelihood(cp,n)   = mvnpdf(y(:,n)',mu,[m_sigma,0;0,m_sigma]);    
                    wp(cp,n)     = wp(cp,n-1)*likelihood(cp,n); % ASSUMES IMPORTANCE DENSITY = PRIOR
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
        xPHatMean(3,n) = sum(wp(:,n).*xp(:,3,n)); % wrapToPi(sum(wp(:,n).*wrapTo2Pi(xp(:,3,n))));
        xPHatMean(4,n) = sum(wp(:,n).*xp(:,4,n));
        xPHatMean(5,n) = sum(wp(:,n).*xp(:,5,n));
    %     xPHatMean = xPHatMean';
    end
    
    t_end = round((cputime - t_start)/60,3);
    fprintf('doin the things a particle can...%.2f\n',t_end);
    %% ERROR
    
    for r = 1:nSamples
%         error(r) = mean((x(1,r)-xPHatMean(1,r))^2);
          error_dist(r) = mean((x(1,1:r)-xPHatMean(1,1:r)).^2 + (x(2,1:r)-xPHatMean(2,1:r)).^2); 
          error_t(r)    = mean((x(3,1:r)-xPHatMean(3,1:r)).^2);
          error_v(r)    = mean((x(4,1:r)-xPHatMean(4,1:r)).^2);
          error_y(r) = mean((x(1,1:r)-y(1,1:r)).^2 + (x(2,1:r)-y(2,1:r)).^2);

    end

    %% PLOT PARTICLES
    
    figure;
    LonP(:,:) = xp(:,1,:);
    plot(LonP','g');
    hold on; 
    plot(y(1,:),'r-.');
    plot(x(1,:),'b-','marker','.');
    plot(xPHatMean(1,:),'k-','marker','.');
    title('Longitude');
    legend('Measured','True','Estimated');

    figure;
    LatP(:,:) = xp(:,2,:);
    plot(LatP','g');
    hold on; plot(x(2,:),'b-','marker','.');
    plot(xPHatMean(2,:),'k-','marker','.');
    title('Latitude');

    figure;
    tP(:,:) = xp(:,3,:);
    plot(tP','g');
    hold on; plot(x(3,:),'b-','marker','.');
    plot(xPHatMean(3,:),'k-','marker','.');
    title('Theta');

    figure;
    vP(:,:) = xp(:,4,:);
    plot(vP','g');
    hold on; plot(x(4,:),'b-','marker','.');
    plot(xPHatMean(4,:),'k-','marker','.');
    title('Velocity');

    figure;
    mP(:,:) = xp(:,5,:);
    plot(mP','g')
    hold on; plot(x(5,:),'b-','marker','.');
    plot(xPHatMean(5,:),'k-','marker','.');
    title('Average Velocity');

    %% PLOT WEIGHTS
    figure;
    subplot(3,1,1);
    plot(0:nSamples-1,std(wp),'r');
    ylabel('$\sigma_w$','Interpreter','LaTeX','FontSize',16);
    xlim([0 nSamples-1]);
    %ylim([0.1^(nSamples) 1]);
    title('weights','Interpreter','LaTeX','FontSize',16);
    box off;

    subplot(3,1,2);
    semilogy(0:nSamples-1,wp,'b');
    ylabel('$w_n^i$','Interpreter','LaTeX','FontSize',16);
    xlim([0 nSamples-1]);
    ylim([0.1^(nSamples) 1]);
    box off;

    subplot(3,1,3);
    plot(0:nSamples-1,wp,'b');
    xlabel('$n$','Interpreter','LaTeX','FontSize',16);
    ylabel('$w_n^i$','Interpreter','LaTeX','FontSize',16);
    xlim([0 nSamples-1]);
    ylim([0 1]);
    box off;

    %% PLOT TRUE VS MEASURED VS MEASURED PATH

    figure;
    plot(y(1,:),y(2,:),'r-o',x(1,:),x(2,:),'b-*',xPHatMean(1,:),xPHatMean(2,:),'g-x');
    legend('Measured','True','Estimated')
%% PLOT POSTERIOR
       
%     PosteriorCalcsPlots(xp,wp,y,y_hat,nRegions,kernelWidth)
    PosteriorCalcsPlots(LonP,wp,y(1,:),xPHatMean(1,:),nRegions,kernelWidth);
    hold on;
    pt = plot(0:99,x(1,:),'b','marker','.');
    hold off;
    xlim([30 50]);
    legend('Measured','Estimated')
    ylabel('Horizontal Position','Interpreter','latex','FontSize',16)
    
    PosteriorCalcsPlots(LatP,wp,y(2,:),xPHatMean(2,:),nRegions,kernelWidth);
    hold on;
    plot(0:99,x(2,:),'b','MarkerSize',10,'marker','.');
    hold off;
    xlim([30 50]);
    legend('Measured','Estimated')
    ylabel('Vertical Position','Interpreter','latex','FontSize',16)
    
    PosteriorCalcsPlots(tP,wp,x(3,:),xPHatMean(3,:),nRegions,kernelWidth);
    hold on;
    plot(0:99,x(3,:),'b','MarkerSize',10,'marker','.');
    hold off;  
    xlim([30 50]);
    legend('True','Estimated')
    ylabel('Heading','Interpreter','latex','FontSize',16)

    PosteriorCalcsPlots(vP,wp,x(4,:),xPHatMean(4,:),nRegions,kernelWidth);
    hold on;
    plot(0:99,x(4,:),'b','MarkerSize',10,'marker','.');
    hold off;
    xlim([30 50]);
    legend('True','Estimated')
    ylabel('Velocity','Interpreter','latex','FontSize',16)

    PosteriorCalcsPlots(mP,wp,x(5,:),xPHatMean(5,:),nRegions,kernelWidth);
    xlim([30 50]);
    legend('True','Estimated')
    ylabel('Mean Velocity','Interpreter','latex','FontSize',16)

    %% STATE SPACE MODEL

    function [xNext] = SystemModel(x,parameters,A,n)
        a = parameters(1);
        b = binornd(1,0.9);
        x1 = x(1)+x(4)*cos(x(3));
        x2 = x(2)+x(4)*sin(x(3));
        x3 = A(x(1),x(2)); % x(3);
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
% end