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
% Ph.bt -- external input
% Ph.Q  -- state noise covarince
% Ph.x0 -- initial state mean
% Ph.Q0 -- initial state covariance
% Ph.C  -- observation matrix
% Ph.d  -- observation constant
% Ph.R  -- observation covariance
% Ph.D  -- observation history effect matrix
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

function [err, y_est] = loo (y, Ph) % [err, y_norm, y_est] = loo (y, Ph)

    [n, T] = size(y);

    x0     = Ph.x0;
    Q0     = Ph.Q0;
    A      = Ph.A;
    Q      = Ph.Q;
    C      = Ph.C;
    R      = Ph.R;
    
    if isfield (Ph, 'd')
        d      = Ph.d;
    else
        d      = zeros(n,1);
    end
    
    m      = length(x0);
    
    if isfield(Ph, 'D')
        D       = Ph.D;
        T       = T - 1;
        yt      = y(:,2:(T+1));
        ys      = y(:,1:T);
        yt      = yt - D * ys;
    else
        yt      = y;
    end
    
    if isfield(Ph, 'bt')
        bt      = Ph.bt;
    else
        bt      = zeros(m,T-1);
    end
    
    err         = 0;
    % y_norm      = 0;
    y_est       = nan(n,T);
    
    for n_neuron   = 1:n
        C_minus_n             = C;
        R_minus_n             = R;
        d_minus_n             = d;
        C_minus_n(n_neuron,:) = [];
        R_minus_n(n_neuron,:) = [];
        d_minus_n(n_neuron,:) = [];
        C_n                   = C(n_neuron,:);
        d_n                   = d(n_neuron,:);
        y_minus_n             = yt;
        y_minus_n(n_neuron,:) = [];
        y_n                   = yt(n_neuron,:);
        Est                   = kalman_forward_backward(y_minus_n, ...
                                A, Q, C_minus_n, R_minus_n, x0, Q0, ...
                                'bt',bt,'d',d_minus_n);
        X_t_minus_n           = Est.X_t;
        Est_y_n               = C_n * X_t_minus_n + d_n;
        err                   = err + norm(Est_y_n - y_n);    
        % y_norm                = y_norm + norm(Est_y_n);
        y_est(n_neuron,:)     = Est_y_n;
    end
    
    if isfield(Ph, 'D')
        ys_est                = ys;
        ys_est(:,2:T)         = y_est(:,1:T-1);
        y_est                 = y_est + D*ys_est;
    end