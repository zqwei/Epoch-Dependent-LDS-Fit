% 1.
% 
% Test of LDS using predefined neural spiking data from a Poission
% distrubtion.
% 
% y      ~ Poiss(exp(z))
% z(t)   = Cx(t)
% x(t+1) = Ax(t) + w(t)
% Estimated parameters include: x0, Q0, A, Q, C
%
%
%  Copyright (C) Ziqiang Wei
%  weiz@janelia.hhmi.org
%  07/22/2014
%
%

clear all %#ok<*CLALL>
addpath('../Release_LDSI_v3/');

% %%
% % 1.
% % 
% This is an example for stable dynamical system ,where |A-1|<0
% Since it is hard to intepret P matrix (Uniqueness of the solution is not
% granteed in this condition) when the dimension of Q is larger than one, 
% the following code is only to consider the case where xDim = 1
% One can futher play with this code for higher dimensionality.
%
% Test with the same y Observation but different possion generation process
%
%
%
rng('Shuffle');

Arot     = 0.1;
Aspec    = 0.99;
Arand    = 0.03;
Q0max    = 0.3;
Rmin     = 0.1;
Rmax     = 10.1;

xDim     = 1;
yDim     = 100;
T        = 80;
nTrial   = 2500;
A        = eye(xDim)+Arand*randn(xDim);
A        = A./max(abs(eig(A)))*Aspec;
MAS      = randn(xDim); MAS = (MAS-MAS')/2;LDS.A  = expm(Arot.*(MAS))*A;
LDS.Q    = eye(xDim);
LDS.V0   = 0.1*eye(xDim);
LDS.x0   = randn(xDim,1)/3;
LDS.C    = rand(yDim,xDim)+0.5./sqrt(3*xDim);
LDS.R    = 0;

%% generate outputs of LDS
% no difference across trials
[X,Z]    = SimulateLDSIdenticalTrials(LDS,T,nTrial);
Yratio   = 1;
Y        = poissrnd(exp(Z)*Yratio);
% Y (Y==0) = 0.001;


Ph       = lds(sqrt(Y), xDim,'mean_type','no_mean','tol',1e-5); 
[err, y_est,~]   = loo(sqrt(Y), Ph);
% 'stage mean' or 'no mean' does no matter for this fit
% P        = sqrt(Ph.Q);

% Result
% disp(['Difference of Latent system (per Chanel, per Time Point, per Trial):   ',...
%     num2str(sqrt(norm(abs(Ph.Xk_t(:))/P - abs(X(:)))^2/norm(X(:))^2/xDim/T/nTrial))]);
1 - err
figure; plot_n_trial = 4;
for n_plot = 1: plot_n_trial
    subplot(1, 4, n_plot)
    hold on;
    % black line: real dynamic in latent space
    % red line  : fitted dynamic in latent space
    plot(sqrt(squeeze(Y(5,:,n_plot))),'--','color',[0.5 0.5 0.5], 'linewid', 1)
    plot(sqrt(exp(squeeze(Z(5,:,n_plot)))*Yratio),'-k', 'linewid', 1)
    plot(squeeze(y_est(5,:,n_plot)),'--r')
    hold off     
    box off
    xlabel('Time')
    ylabel('Sqrt. root simulated firing rate')
end
setPrint(4*8, 6, 'Plots/Test_Poisson_High_firing_rate_Idential_Trial')

Yratio   = 0.5;
Y        = poissrnd(exp(Z)*Yratio);
Ph       = lds(sqrt(Y), xDim,'mean_type','no_mean','tol',1e-5); 
[err, y_est,~]   = loo(sqrt(Y), Ph);
1 - err
figure; plot_n_trial = 4;
for n_plot = 1: plot_n_trial
    subplot(1, 4, n_plot)
    hold on;
    plot(sqrt(squeeze(Y(5,:,n_plot))),'--','color',[0.5 0.5 0.5], 'linewid', 1)
    plot(sqrt(exp(squeeze(Z(5,:,n_plot)))*Yratio),'-k', 'linewid', 1)
    plot(squeeze(y_est(5,:,n_plot)),'--r')
    hold off
    box off
    xlabel('Time')
    ylabel('Sqrt. root simulated firing rate')
end
setPrint(4*8, 6, 'Plots/Test_Poisson_Medium_firing_rate_Idential_Trial')

Yratio   = 0.1;
Y        = poissrnd(exp(Z)*Yratio);
Ph       = lds(sqrt(Y), xDim,'mean_type','no_mean','tol',1e-5); 
[err, y_est,~]   = loo(sqrt(Y), Ph);
1 - err
figure; plot_n_trial = 4;
for n_plot = 1: plot_n_trial
    subplot(1, 4, n_plot)
    hold on;
    plot(sqrt(squeeze(Y(5,:,n_plot))),'--','color',[0.5 0.5 0.5], 'linewid', 1)
    plot(sqrt(exp(squeeze(Z(5,:,n_plot)))*Yratio),'-k', 'linewid', 1)
    plot(squeeze(y_est(5,:,n_plot)),'--r')
    hold off
    box off
    xlabel('Time')
    ylabel('Sqrt. root simulated firing rate')
end
setPrint(4*8, 6, 'Plots/Test_Poisson_Low_firing_rate_Idential_Trial')


%% generate outputs of LDS
% no difference across trials
[X,Z]    = SimulateLDS(LDS,T,nTrial);
Z        = Z/max(Z(:))*4.3;
Y        = poissrnd(exp(Z)*10);
Ph       = lds(sqrt(Y), xDim,'mean_type','no_mean','tol',1e-5); 
[~, y_est,~]   = loo(sqrt(Y), Ph);
% 'stage mean' or 'no mean' does no matter for this fit
% P        = sqrt(Ph.Q);

% Result
% disp(['Difference of Latent system (per Chanel, per Time Point, per Trial):   ',...
%     num2str(sqrt(norm(abs(Ph.Xk_t(:))/P - abs(X(:)))^2/norm(X(:))^2/xDim/T/nTrial))]);

figure; plot_n_trial = 4;
for n_plot = 1: plot_n_trial
    subplot(1, 4, n_plot)
    hold on;
    % black line: real dynamic in latent space
    % red line  : fitted dynamic in latent space
    plot(sqrt(squeeze(Y(5,:,n_plot))),'--','color',[0.5 0.5 0.5], 'linewid', 1)
    plot(sqrt(exp(squeeze(Z(5,:,n_plot)))*10),'--k', 'linewid', 1)
    plot(squeeze(y_est(5,:,n_plot)),'--r')
    hold off     
    box off
    xlabel('Time')
    ylabel('Sqrt. root simulated firing rate')
end
setPrint(4*8, 6, 'Plots/Test_Poisson_High_firing_rate_Random_Trial')

Y        = poissrnd(exp(Z));
Ph       = lds(sqrt(Y), xDim,'mean_type','no_mean','tol',1e-5); 
[~, y_est,~]   = loo(sqrt(Y), Ph);
figure; plot_n_trial = 4;
for n_plot = 1: plot_n_trial
    subplot(1, 4, n_plot)
    hold on;
    plot(sqrt(squeeze(Y(5,:,n_plot))),'--','color',[0.5 0.5 0.5], 'linewid', 1)
    plot(sqrt(exp(squeeze(Z(5,:,n_plot)))),'--k', 'linewid', 1)
    plot(squeeze(y_est(5,:,n_plot)),'--r')
    hold off
    box off
    xlabel('Time')
    ylabel('Sqrt. root simulated firing rate')
end
setPrint(4*8, 6, 'Plots/Test_Poisson_Medium_firing_rate_Random_Trial')

Y        = poissrnd(exp(Z)*0.1);
Ph       = lds(sqrt(Y), xDim,'mean_type','no_mean','tol',1e-5); 
[~, y_est,~]   = loo(sqrt(Y), Ph);
figure; plot_n_trial = 4;
for n_plot = 1: plot_n_trial
    subplot(1, 4, n_plot)
    hold on;
    plot(sqrt(squeeze(Y(5,:,n_plot))),'--','color',[0.5 0.5 0.5], 'linewid', 1)
    plot(sqrt(exp(squeeze(Z(5,:,n_plot)))*0.1),'--k', 'linewid', 1)
    plot(squeeze(y_est(5,:,n_plot)),'--r')
    hold off
    box off
    xlabel('Time')
    ylabel('Sqrt. root simulated firing rate')
end
setPrint(4*8, 6, 'Plots/Test_Poisson_Low_firing_rate_Random_Trial')



%%
% % % This is an example for stable dynamical system ,where |A-1|<0
% % % Since it is hard to intepret P matrix (Uniqueness of the solution is not
% % % granteed in this condition) when the dimension of Q is larger than one, 
% % % the following code is only to consider the case where xDim = 1
% % % One can futher play with this code for higher dimensionality.
% % %
% % rng('Shuffle');
% % 
% % Arot     = 0.1;
% % Aspec    = 0.99;
% % Arand    = 0.03;
% % Q0max    = 0.3;
% % Rmin     = 0.1;
% % Rmax     = 10.1;
% % 
% % xDim     = 1;
% % yDim     = 100;
% % T        = 80;
% % nTrial   = 2500;
% % A        = eye(xDim)+Arand*randn(xDim);
% % A        = A./max(abs(eig(A)))*Aspec;
% % MAS      = randn(xDim); MAS = (MAS-MAS')/2;LDS.A  = expm(Arot.*(MAS))*A;
% % LDS.Q    = eye(xDim);
% % LDS.V0   = 0.1*eye(xDim);
% % LDS.x0   = randn(xDim,1)/3;
% % LDS.C    = randn(yDim,xDim)./sqrt(3*xDim);
% % LDS.R    = 0;
% % 
% % % generate outputs of LDS
% % [X,Z]    = SimulateLDS(LDS,T,nTrial);
% % Y        = poissrnd(exp(Z));
% % % Y (Y==0) = 0.001;
% % 
% % 
% % Ph       = lds(log10(Y), xDim,'mean_type','no_mean','tol',1e-5); 
% % % 'stage mean' or 'no mean' does no matter for this fit
% % P        = sqrt(Ph.Q);
% % 
% % % Result
% % disp(['Difference of Latent system (per Chanel, per Time Point, per Trial):   ',...
% %     num2str(sqrt(norm(abs(Ph.Xk_t(:))/P - abs(X(:)))^2/norm(X(:))^2/xDim/T/nTrial))]);
% % 
% % figure; plot_n_trial = 16;
% % for n_plot = 1: plot_n_trial
% %     subplot(4, 4, n_plot)
% %     hold on;
% %     % black line: real dynamic in latent space
% %     % red line  : fitted dynamic in latent space
% %     plot(abs(squeeze(X(:,:,n_plot)))*P,'-k')
% %     plot(abs(squeeze(Ph.Xk_t(:,:,n_plot))),'--r')
% %     hold off        
% % end
% % 
% % disp(['Difference of parameter A (nomalized by norm(A)):   ', num2str(norm(abs(LDS.A)-abs(Ph.A))/norm(LDS.A))]);
% % figure; hold on; plot(Ph.C, LDS.C,'ok');plot([-1 1],[-1/P 1/P],'-r'); plot([-1 1],[1/P -1/P],'-r'); hold off
% % xlabel('Estimated matrix of C');
% % ylabel('True matrix of C')

% %%
% % 2.
% % 
% % Test of LDS using predefined neural spiking data from a rounded Lognormal
% % distrubtion.
% % 
% % y      = round(exp(z))
% % z(t)   = Cx(t) + v(t)
% % x(t+1) = Ax(t) + w(t)
% % Estimated parameters include: x0, Q0, A, Q, C
% %
% %
% %  Copyright (C) Ziqiang Wei
% %  weiz@janelia.hhmi.org
% %  07/22/2014
% %
% %
% 
% clear all
% addpath('../Release_LDSI_v3/');
% 
% % This is an example for stable dynamical system ,where |A-1|<0
% % Since it is hard to intepret P matrix (Uniqueness of the solution is not
% % granteed in this condition) when the dimension of Q is larger than one, 
% % the following code is only to consider the case where xDim = 1
% % One can futher play with this code for higher dimensionality.
% %
% rng('Shuffle');
% 
% Arot     = 0.1;
% Aspec    = 0.99;
% Arand    = 0.03;
% Q0max    = 0.3;
% Rmin     = 0.1;
% Rmax     = 10.1;
% 
% xDim     = 1;
% yDim     = 100;
% T        = 80;
% nTrial   = 2500;
% A        = eye(xDim)+Arand*randn(xDim);
% A        = A./max(abs(eig(A)))*Aspec;
% MAS      = randn(xDim); MAS = (MAS-MAS')/2;
% LDS.A    = expm(Arot.*(MAS))*A;
% LDS.A    = 0.7;
% LDS.Q    = eye(xDim);
% LDS.V0   = 0.1*eye(xDim);
% LDS.x0   = randn(xDim,1)/3;
% LDS.C    = randn(yDim,xDim)./sqrt(3*xDim);
% LDS.R    = diag(rand(yDim,1)*Rmax+Rmin);
% 
% % generate outputs of LDS
% [X,Z]    = SimulateLDS(LDS,T,nTrial);
% Y        = round(exp(Z));
% Y (Y==0) = 0.01; % This value should be adjusted, it is not a fixed one.
% 
% 
% Ph       = lds(log10(Y), xDim,'mean_type','no_mean','tol',1e-5); 
% % 'stage mean' or 'no mean' does no matter for this fit
% P        = sqrt(Ph.Q);
% 
% % Result
% disp(['Difference of Latent system (per Chanel, per Time Point, per Trial):   ',...
%     num2str(sqrt(norm(abs(Ph.Xk_t(:))/P - abs(X(:)))^2/norm(X(:))^2/xDim/T/nTrial))]);
% 
% figure; plot_n_trial = 16;
% for n_plot = 1: plot_n_trial
%     subplot(4, 4, n_plot)
%     hold on;
%     % black line: real dynamic in latent space
%     % red line  : fitted dynamic in latent space
%     plot(abs(squeeze(X(:,:,n_plot)))*P,'-k')
%     plot(abs(squeeze(Ph.Xk_t(:,:,n_plot))),'--r')
%     hold off   
%     box off
%     xlabel('Time')
%     ylabel('Simulated data')
% end
% title('Single-unit LDS model fit')


% disp(['Difference of parameter A (nomalized by norm(A)):   ', num2str(norm(abs(LDS.A)-abs(Ph.A))/norm(LDS.A))]);
% figure; hold on; plot(Ph.C, LDS.C,'ok');plot([-1 1],[-1/P 1/P],'-r'); plot([-1 1],[1/P -1/P],'-r'); hold off
% xlabel('Estimated matrix of C');
% ylabel('True matrix of C')



