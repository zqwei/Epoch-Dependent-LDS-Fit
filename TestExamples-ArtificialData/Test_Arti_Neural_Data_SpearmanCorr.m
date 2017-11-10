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
addpath('../Code/');


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

% no difference across trials
[X,Z]    = SimulateLDS(LDS,T,nTrial);
Z        = Z/max(Z(:))*4.3;
Y        = poissrnd(exp(Z)*10);
Ph       = lds(sqrt(Y), xDim,'mean_type','no_mean','tol',1e-5); 
[~, y_est,~]   = loo(sqrt(Y), Ph);
disp(['High firing rate; Random Trial: ' num2str(corr(sqrt(Y(:)), y_est(:), 'type', 'Spearman'))])
disp(['High firing rate; Random Trial; expectation: ' num2str(corr(sqrt(exp(Z(:))*10), y_est(:), 'type', 'Spearman'))])

Y        = poissrnd(exp(Z));
Ph       = lds(sqrt(Y), xDim,'mean_type','no_mean','tol',1e-5); 
[~, y_est,~]   = loo(sqrt(Y), Ph);
disp(['Medium firing rate; Random Trial: ' num2str(corr(sqrt(Y(:)), y_est(:), 'type', 'Spearman'))])
disp(['Medium firing rate; Random Trial; expectation: ' num2str(corr(sqrt(exp(Z(:))), y_est(:), 'type', 'Spearman'))])

Y        = poissrnd(exp(Z)*0.1);
Ph       = lds(sqrt(Y), xDim,'mean_type','no_mean','tol',1e-5); 
[~, y_est,~]   = loo(sqrt(Y), Ph);
disp(['Low firing rate; Random Trial; raw data: ' num2str(corr(sqrt(Y(:)), y_est(:), 'type', 'Spearman'))])
disp(['Low firing rate; Random Trial; expectation: ' num2str(corr(sqrt(exp(Z(:))*0.1), y_est(:), 'type', 'Spearman'))])




