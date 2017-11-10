%
% [A, Q, bt] = M_step_state (X_t, Q_xx_ts, Q_xx_ss, Q_xx_uu, K, is_bt, is_diag)
%
% Parameter estimation for linear dynamic system (M step)
%
% Inputs:
% 
% X_t      --  sum_k (E(x_t))
% Q_xx_ts  --  sum_k sum_t (E(x_t * x_s'))
% Q_xx_ss  --  sum_k sum_t (E(x_s * x_s'))
% Q_xx_uu  --  sum_k sum_t (E(x_t * x_t'))
% K        --  number of trials
% is_bt    --  true: if model includes bt
%              default: false
% is_diag  --  true: if matrix Q is diagonal (white noise)
%              default: true
% sum_t    --  from t = 2 to T
%
% Output:
% 
% x0    -- initial state mean
% Q0    -- initial state covariance
%
%
% Model:
%
%             y(k,t) = Ph.C * x(k,t) + Ph.D * y(k,s) + Ph.d(k) + v(k,t)
%             x(k,t) = Ph.A * x(k,s) + Ph.bt(k,s) + w(k,s)
%             s      = t - 1
%        where
%             v ~ N(0,R)
%             w ~ N(0,Q)
%             x(k,1) ~ N(pi,Q0) (for any k)
%
% Ver: 2.0
%
% @ 2014 Ziqiang Wei
% weiz@janelia.hhmi.org
% 
% 
function [A, Q, is_success] = M_step_state (Q_xx_ts, Q_xx_ss, Q_xx_uu, TK)
    
    minVar     = 1e-8;    
    is_success = true;
    invQ_xx_ss = invSymFast(Q_xx_ss);
    A          = Q_xx_ts * invQ_xx_ss;
    Q          = 1/TK * (Q_xx_uu - A*Q_xx_ts');
    
    Q          = diag(Q);
    
    if any(Q)  <=0
        is_success = false;
        display('Using private minimum of Q!');
        Q_max      = max(Q);
        Q(Q<=0)    = max(Q_max * minVar, minVar);
    end
    Q          = diag(Q);
