clear all
addpath('../Code/');

% This is an example for stable dynamical system ,where |A-1|<0
% Since it is hard to intepret P matrix (Uniqueness of the solution is not
% granteed in this condition) when the dimension of Q is larger than one, 
% the following code is only to consider the case where xDim = 1
% One can futher play with this code for higher dimensionality.
%
rng('Shuffle');

Arot     = 0.1;
Aspec    = 0.99;
Arand    = 0.03;
Q0max    = 0.3;
Rmin     = 0.1;
Rmax     = 10.1;

xDim     = 3;
yDim     = 10;
T        = 80;
nTrial   = 2500;
A        = eye(xDim)+Arand*randn(xDim);
A        = A./max(abs(eig(A)))*Aspec;
MAS      = randn(xDim); MAS = (MAS-MAS')/2;
LDS.A    = expm(Arot.*(MAS))*A;
LDS.Q    = eye(xDim);
LDS.V0   = 0.1*eye(xDim);
LDS.x0   = randn(xDim,1)/3;
LDS.C    = randn(yDim,xDim)./sqrt(3*xDim);
LDS.R    = diag(rand(yDim,1)*Rmax+Rmin);

% generate outputs of LDS
[X,Y, Y_est]    = SimulateLDS(LDS,T,nTrial);

Per_exp  = zeros(1,4);
BIC      = zeros(1,4);
AIC      = zeros(1,4);

for xDim           = 1:5;
    Ph             = lds(Y, xDim,'mean_type','no_mean','tol',1e-5); 
    % 'stage mean' or 'no mean' does no matter for this fit
    P              = sqrt(Ph.Q);
    Y_all          = reshape(Y, yDim, []);
    R0             = cov(Y_all');
    Per_exp(xDim)  = 1-sum(diag(Ph.R))/sum(diag(R0));
    KC             = yDim * xDim;
    KR             = yDim;
    KA             = xDim*xDim;
    KQ             = xDim;
    BIC(xDim)      = -2*nanmax(Ph.LL) + (KC+KR+KA+2*KQ+xDim)*log(yDim*T*nTrial);
    AIC(xDim)      = -2*nanmax(Ph.LL) + 2*(KC+KR+KA+2*KQ+xDim);
end

figure;
subplot(1, 3, 1)
plot(Per_exp,'-ok','linewid', 2); ylabel('% Explained Variance')
box off; xlabel('# of Latent Dimensions')
subplot(1, 3, 2)
plot(BIC/max(BIC),'-ok','linewid', 2); ylabel('Normalized BIC')
box off; xlabel('# of Latent Dimensions')
subplot(1, 3, 3)
plot(AIC/max(AIC),'-ok','linewid', 2); ylabel('Normalized AIC')
box off; xlabel('# of Latent Dimensions')
setPrint(3*8, 6, 'Plots/Test_multi_latent_unit_performance', 'png')

Ph             = lds(Y, 3,'mean_type','no_mean','tol',1e-5);
[~, y_est,~]   = loo(Y, Ph);
figure; plot_n_trial = 9;
for n_plot = 1: plot_n_trial
    subplot(3, 3, n_plot)
    hold on;
    % black line: real dynamic in latent space
    % red line  : fitted dynamic in latent space
    plot(squeeze(Y(ceil(n_plot/3),:,n_plot)),'--','color',[0.5 0.5 0.5])
    plot(squeeze(Y_est(ceil(n_plot/3),:,n_plot)),'-k','linewid',1)
    plot(squeeze(y_est(ceil(n_plot/3),:,n_plot)),'--r','linewid',1)
    hold off  
    box off
    xlabel('Time')
    ylabel('Simulated observed data')
end

suptitle('Three latent units')
setPrint(3*8, 3*6, 'Plots/Test_multi_latent_unit', 'png')

% err_all            = sum(Y(:).^2);
% 
% for xDim           = 1:4;
%     Ph             = lds(Y, xDim,'mean_type','no_mean','tol',1e-5); 
%     [~, y_est,~]   = loo(Y, Ph);
%     Per_exp(xDim)  = 1-sum((y_est(:)-Y(:)).^2)/err_all;
% end
% 
% est_exp            = 1-sum((Y_est(:)-Y(:)).^2)/err_all;

% err = cross_valid_ldsi(Y, 4, 4,[], 'mean_type','no_mean','tol',1e-5);
% disp(Per_exp' - (1- err/sum(diag(R0))/T/nTrial*4))
% These two variables should be nearly the same.