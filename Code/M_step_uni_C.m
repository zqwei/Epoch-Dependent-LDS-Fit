%
% [C, R, D] = ...
%    M_step_obs (Q_yy_tt, Q_yx_tt, Q_xx_tt, TK, is_diag, is_D, Q_yy_ts, Q_xy_ts, invQ_yy_ss)
%
% Parameter estimation for linear dynamic system (M step)
%
% Inputs:
% 
% X_t      --  sum_k (E(x_t))
% Q_yy_tt  --  sum_k sum_t (E(y_t * y_t'))
% Q_yx_tt, Q_xx_tt, Q_yy_ts, Q_xy_ts, invQ_yy_ss follow the similar definition 
% TK       --  multiplication of number of trials and length of time series
% is_D     --  true: if model includes D
%              default: false
% is_diag  --  true: if matrix R is diagonal (white noise)
%              default: true
% sum_t    --  from t = 1 to T
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
% Ver: 1.0
%
% @ 2014 Ziqiang Wei
% weiz@janelia.hhmi.org
% 
% 
function C = M_step_uni_C (Q_yx_tt, Q_xx_tt)
    invQ_xx_tt = invSymFast(Q_xx_tt);
    C          = Q_yx_tt * invQ_xx_tt;
