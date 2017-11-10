% Test of LDS
% Following - 2005 Sen Cheng
% This test code shows
%

%%
% 1. whether we can uncover the true low dimensional subspace in the limit
%    of large data for single stage.
%    - It is not quite reasonable to compare matrix A and Q directly with
%      the true parameters, since any invertible matrix P can lead to the
%      indentical observation as long as:
%      xt -> Pxt
%      yt -> yt
%      A  -> P*A*invP (if Q is one-dimension, A -> A)
%      B  -> P*B
%      Q  -> P*Q*Pt   (if Q is one-dimension, Q -> Q*p^2)
%      C  -> C*invP   (if Q is one-dimension, C -> C/p)
%      D  -> D
%      R  -> R
%      (See Cheng & Sabes, 2006, Neural Comput.)
%    - At the same time, y_est is never equal to y_t, due to the existence
%      of R.
%    - One solution is fixed the dim of Q as 1 or the value of Q in updatae, 
%      say eye(xdim) for instance, and see whether x_est converages to x_t. 
%      (See Pfau, Pnevmatikakis, Paninski, 2013, NIPS)
%
%
%  Copyright (C) Ziqiang Wei
%  weiz@janelia.hhmi.org
%  07/22/2014
%
%

% clear all
% addpath('../Release_LDSI_v3/');
% 
% % 1.1
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
% MAS      = randn(xDim); MAS = (MAS-MAS')/2;LDS.A  = expm(Arot.*(MAS))*A;
% LDS.Q    = eye(xDim);
% LDS.V0   = 0.1*eye(xDim);
% LDS.x0   = randn(xDim,1)/3;
% LDS.C    = randn(yDim,xDim)./sqrt(3*xDim);
% LDS.R    = diag(rand(yDim,1)*Rmax+Rmin);
% 
% % generate outputs of LDS
% [X,Y]  = SimulateLDS(LDS,T,nTrial);
% 
% Ph     = lds(Y, xDim,'mean_type','no_mean','tol',1e-5); 
% % 'stage mean' or 'no mean' does no matter for this fit
% P      = sqrt(Ph.Q);
% 
% % Here I also put a SSID implementation of the code
% % System from SSID
% %    x(t+1)      = Ax(t) + Ke(t)
% %    y(t)        = Cx(t) + e(t)
% %    where R : = cov(e(t)) & Q := K*R*K'
% % This might result from lacking of estimation of x0 in the system
% % SSID is thus a good way to initilize parameter like A and C
% %                         to estimate a minimum required latent dimension
% % 
% % I don't feel it is quite useful for parameter Q, R and Q0
% [A, C, K, R] = stochastic(Y, 1, 2);
% 
% 
% % Result
% disp(['Difference of Latent system (per Chanel, per Time Point, per Trial):   ',...
%     num2str(sqrt(norm(abs(Ph.Xk_t(:))/P - abs(X(:)))^2/norm(X(:))^2/xDim/T/nTrial))]);
% disp(['Difference of parameter A (nomalized by norm(A)):   ', num2str(norm(abs(LDS.A)-abs(Ph.A))/norm(LDS.A))]);
% disp(['Difference of parameter C (nomalized by norm(C)):   ', num2str(norm(abs(LDS.C)-abs(Ph.C)*P)/norm(LDS.C))]);
% disp(['Difference of parameter R (nomalized by norm(R)):   ', num2str(norm(abs(LDS.R)-abs(Ph.R))/norm(LDS.R))]);
% 
% % 1.2
% % Following loop is to redo the same computation above but with multiple
% % times and get some statistical result for the reliablility of the fit.
% totFit   = 100;
% err_X    = nan(totFit,1);
% err_A    = nan(totFit,1);
% err_C    = nan(totFit,1);
% err_R    = nan(totFit,1);
% nRmax    = nan(totFit,1);
% 
% for nFit = 1:totFit
%     Arot     = 0.1;
%     Aspec    = 0.99;
%     Arand    = 0.03;
%     Q0max    = 0.3;
%     Rmin     = 0.1;
%     Rmax     = rand()*100+0.1;
%     nRmax(nFit) = Rmax;
% 
%     xDim     = 1;
%     yDim     = 100;
%     T        = 80;
%     nTrial   = 2500;
%     A        = eye(xDim)+Arand*randn(xDim);
%     A        = A./max(abs(eig(A)))*Aspec;
%     MAS      = randn(xDim); MAS = (MAS-MAS')/2;LDS.A  = expm(Arot.*(MAS))*A;
%     LDS.Q    = eye(xDim);
%     LDS.V0   = 0.1*eye(xDim);
%     LDS.x0   = randn(xDim,1)/3;
%     LDS.C    = randn(yDim,xDim)./sqrt(3*xDim);
%     LDS.R    = diag(rand(yDim,1)*Rmax+Rmin);
% 
%     % generate outputs of LDS
%     [X,Y]  = SimulateLDS(LDS,T,nTrial);
% 
%     Ph     = lds(Y, xDim);
%     P      = sqrt(Ph.Q);
% 
%     % Result
%     err_X(nFit) = (norm(abs(Ph.Xk_t(:))/P - abs(X(:)))^2)/(norm(X(:))^2)/xDim/T/nTrial;
%     err_A(nFit) = norm(abs(LDS.A)-abs(Ph.A))/norm(LDS.A);
%     err_C(nFit) = norm(abs(LDS.C)-abs(Ph.C)*P)/norm(LDS.C);
%     err_R(nFit) = norm(abs(LDS.R)-abs(Ph.R))/norm(LDS.R);
% end
% 
% % The following figure shows that what happens to a different extends of R
% % plot(nRmax,[err_X,err_A,err_C,err_R],'.')
% % The difference in my simulation for each parameter is less than 0.8%


%%
% 2. whether we can uncover the true low dimensional subspace in the limit
%    of large data for multiple stage.
%
%
%  Copyright (C) Ziqiang Wei
%  weiz@janelia.hhmi.org
%  07/22/2014

clear all
addpath('../Code/');

rng('Shuffle');

totStage   = 4; 
Arot       = 0.1;
Aspec      = 0.99;
Arand      = 0.03;
Q0max      = 0.3;
Rmin       = 0.1;
Rmax       = 10.1;

xDim       = 1;
yDim       = 100;
T          = 100;
nTrial     = 2500;
timePoints = [40  60 80];
LDS.A      = zeros(xDim, xDim, totStage);
LDS.Q      = zeros(xDim, xDim, totStage);
LDS.V0     = 0.1*eye(xDim);
LDS.x0     = randn(xDim,1)/3;
LDS.C      = zeros(yDim, xDim, totStage);
LDS.R      = zeros(yDim, yDim, totStage);

for nStage            = 1:totStage
    A                 = eye(xDim)+Arand*randn(xDim);
    A                 = A./max(abs(eig(A)))*Aspec;
    MAS               = randn(xDim); 
    MAS               = (MAS-MAS')/2;
    LDS.A(:,:,nStage) = expm(Arot.*(MAS))*A;
    LDS.Q(:,:,nStage) = eye(xDim);
    LDS.C(:,:,nStage) = randn(yDim,xDim)./sqrt(3*xDim);
    LDS.R(:,:,nStage) = diag(rand(yDim,1)*Rmax+Rmin);
end


% generate outputs of LDS
[X,Y]  = SimulateLDS_MultiStage(LDS, T, nTrial, timePoints);
is_fit = false;
while ~is_fit
    try
        Ph     = lds(Y, xDim, 'timePoint', timePoints, 'mean_type','no_mean');
        is_fit = true;
    catch
        is_fit = false;
    end    
end
P_all        = sqrt(Ph.Q);

timePoint    = [0, timePoints, T];

plot_n_trial = 9;

figure;
for nStage   = 1:totStage
    X_nStage = X(:,timePoint(nStage)+1:timePoint(nStage+1),:);
    X_est    = Ph.Xk_t(:,timePoint(nStage)+1:timePoint(nStage+1),:);
    P        = P_all(:,:,nStage);
    disp(['Stage:    ',num2str(nStage)]);
    disp(['Difference of Latent system (per Chanel, per Time Point, per Trial):   ',...
        num2str(norm(abs(X_est(:))/P - abs(X_nStage(:)))^2/norm(X_nStage(:))^2/xDim/(timePoint(nStage+1)-timePoint(nStage))/nTrial)]);
    disp(['Difference of parameter A (nomalized by norm(A)):   ', num2str(norm(abs(LDS.A(:,:,nStage))-abs(Ph.A(:,:,nStage)))/norm(LDS.A(:,:,nStage)))]);
    disp(['Difference of parameter C (nomalized by norm(C)):   ', num2str(norm(abs(LDS.C(:,:,nStage))-abs(Ph.C(:,:,nStage))*P)/norm(LDS.C(:,:,nStage)))]);
    disp(['Difference of parameter R (nomalized by norm(R)):   ', num2str(norm(abs(LDS.R(:,:,nStage))-abs(Ph.R(:,:,nStage)))/norm(LDS.R(:,:,nStage)))]);
    
    for n_plot = 1: plot_n_trial
        subplot(3, 3, n_plot)
        hold on;
        % black line: real dynamic in latent space
        % red line  : fitted dynamic in latent space
        plot(timePoint(nStage)+1:timePoint(nStage+1),abs(squeeze(X_nStage(:,:,n_plot)))*P,'-k', 'linewid', 1)
        plot(timePoint(nStage)+1:timePoint(nStage+1),abs(squeeze(X_est(:,:,n_plot))),'--r', 'linewid', 1)
        hold off
        box off
        xlabel('Time')
        ylabel('Simulated observed data')
    end
    
end
suptitle('Single latent unit')
setPrint(3*8, 3*6, 'Plots/Test_single_latent_unit', 'png')
