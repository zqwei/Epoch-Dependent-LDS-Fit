%
% Est = kalman_forward_backward(Y, m, lambda, cyc, tol)
%
% Kalman forward backward step (E-step) of parameter estimation for trials
% with the identical Ts
%
% Inputs:
%
% Y   -- n x T x k observation data matrix
% n   -- dimension of observation
% T   -- length of obervations in each trial
%
% Parameters used in the basic LDS system:
% A, Q, C, R, x0, Q0
% 
% optional parameters:
% D (also with ys), d, bt
% when D is provided, ys is also required
%
% Output Est (struct) ->
% 
% Est.X_t        -- estimation of state variable
% Est.Q_xx_tt    -- estimation of variance of state variable
% Est.Q_xx_ts    -- estimation of covariance of state variable
% Est.is_success -- return if E-step is done successfully
% Est.lik        -- likelihood
%
%
% Model:
%
%             y(k,t) = Ph.C * x(k,t) + v(k,t)
%             x(k,t) = Ph.A * x(k,s) + w(k,s)
%             s      = t - 1
%        where
%             v ~ N(0,R)
%             w ~ N(0,Q)
%             x(k,1) ~ N(pi,Q0) (for any k)
%
%
%
% This is main function that handles EM iterations until change of LL < tol
% or cyc step reaches maximum cyc value.
%
% @ 2014 Ziqiang Wei
% weiz@janelia.hhmi.org
% 
% version 2.0

function Est = kalman_forward_backward_xy(y, Ph, timePoint)

    [yDim, T, K] = size(y);
    
    if nargin == 2
        timePoint = [0, T];
    end

    x0           = Ph.x0;
    Q0           = Ph.Q0;
    A            = Ph.A;
    Q            = Ph.Q;
    C            = Ph.C;
    R            = Ph.R;
    d            = Ph.d;
    if size(d,1) ~= yDim
        d        = ones(yDim,1) * d;
    end
    
    if nargin == 2
        nt       = 1;
    else
        nt       = size(C,3);
    end
    
    xDim         = size(A,1);
    
    y_est        = nan(yDim,T,K);
    w            = nan(xDim,T-1,K); % noise parameter (latent network layer)
    
    for nt_now   = 1:nt
        y(:,timePoint(nt_now)+1:timePoint(nt_now+1),:) = ...
            remove_mean(y(:,timePoint(nt_now)+1:timePoint(nt_now+1),:),d(:,nt_now));
    end
    
    Est          = kalman_forward_backward(y, A, Q, C, R, x0, Q0, timePoint);
    
    X_t          = Est.Xk_t;
    
    X_u          = X_t(:, 2:end ,  :);
    X_s          = X_t(:, 1:end-1 ,:);
    timePointx   = timePoint;
    timePointx(end) = timePointx(end) - 1;
    
    for nt_now   = 1:nt
        for nTrial = 1:K
            C_nt                = C(:,:,nt_now);
            A_nt                = A(:,:,nt_now);
            X_nt                = X_t(:,timePoint(nt_now)+1:timePoint(nt_now+1),nTrial);
            y_est(:,timePoint(nt_now)+1:timePoint(nt_now+1),nTrial) = ....
                C_nt * X_nt;
            X_nu                = X_u(:,timePointx(nt_now)+1:timePointx(nt_now+1),nTrial);
            X_ns                = X_s(:,timePointx(nt_now)+1:timePointx(nt_now+1),nTrial);
            
            w(:,timePointx(nt_now)+1:timePointx(nt_now+1),nTrial) = ....
                X_nu - A_nt * X_ns;
        end
    end
    
    Est.v        = y - y_est;
    Est.w        = w;
    