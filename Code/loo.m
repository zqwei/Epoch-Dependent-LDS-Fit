%
% [err, y_norm, y_est] = loo (y, Ph)
% 
% Leave-one-neuron out evaluation of fitted parameters
%
% Inputs:
% 
% y     -- observed data
% Ph    -- parameter set
% Ph.A  -- state transition matrix
% Ph.Q  -- state noise covarince
% Ph.x0 -- initial state mean
% Ph.Q0 -- initial state covariance
% Ph.C  -- observation matrix
% Ph.R  -- observation covariance
%
%
% Output:
% 
% err    -- sum squared estimation error
% y_norm -- sum norm of estimated data
% y_est  -- estimation of y
%
%
% Model:
%
%             y(t) = Ph.C * x(t) + Ph.D * y(s) + Ph.d + v(t)
%             x(t) = Ph.A * x(s) + Ph.bt(s) + w(s)
%             s      = t - 1
%        where
%             v ~ N(0,R)
%             w ~ N(0,Q)
%             x(k,1) ~ N(pi,Q0) (for any k)
%
% leave one neuron out is performed on:
%
%             z(t) = y(t) - D*y(s)
% where
%             z_est(t,i) = Ph.C_i * x_est(t) + Ph.d_i
%             x_est(t)   = E(x|Ph.C_-i, Ph.d_-i, A, bt, R_-i, Q, x0, Q0, z(t, -i))
%
% Ver: 1.0
%
% @ 2014 Ziqiang Wei
% weiz@janelia.hhmi.org
% 
% 

function [err, y_est, rand_y] = loo (y, Ph, timePoint) % Leave one neuron out evaluation
    
    % y --- yDim x T x K

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
    
    
    err          = 0;
    y_est        = nan(yDim,T,K);
    Est_y_n      = nan(1,T,K);
    
    for nt_now   = 1:nt
        y(:,timePoint(nt_now)+1:timePoint(nt_now+1),:) = ...
            remove_mean(y(:,timePoint(nt_now)+1:timePoint(nt_now+1),:),d(:,nt_now));
    end
    
    rand_y       = sum(y(:).^2); % all variance needs to explain
    
    for n_neuron   = 1:yDim
        
        C_minus_n               = C;
        R_minus_n               = R;
        C_minus_n(n_neuron,:,:) = [];
        R_minus_n(n_neuron,:,:) = [];
        R_minus_n(:,n_neuron,:) = [];
        C_n                     = C(n_neuron,:,:);
        y_minus_n               = y;
        y_minus_n(n_neuron,:,:) = [];
        y_n                     = y(n_neuron,:,:);
        
        Est                     = kalman_forward_backward_w_n_lik(y_minus_n, ...
                                    A, Q, C_minus_n, R_minus_n, x0, Q0, ...
                                    timePoint);
        X_t_minus_n             = Est.Xk_t;
        
        for nt_now  = 1:nt
            X_nt                = X_t_minus_n(:,timePoint(nt_now)+1:timePoint(nt_now+1),:);
            X_nt                = reshape(X_nt, xDim, []);
            C_nt                = C_n(:,:,nt_now);
            Est_y_nt            = C_nt * X_nt;
            Est_y_n(:,timePoint(nt_now)+1:timePoint(nt_now+1),:) = reshape(Est_y_nt,1,[],K);
        end
        
        err                     = err + norm(Est_y_n(:) - y_n(:)).^2;
        for nt_now   = 1:nt
            y_est(n_neuron,timePoint(nt_now)+1:timePoint(nt_now+1),:) = ...
                Est_y_n(:,timePoint(nt_now)+1:timePoint(nt_now+1),:) + d(n_neuron,nt_now);
        end
    end
    
    err = err/rand_y;
    