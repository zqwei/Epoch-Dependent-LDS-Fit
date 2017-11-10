%
% [x0, Q0] = M_step_init (x1, Q1, K)
%
% Parameter estimation for linear dynamic system (M step)
%
% Inputs:
% 
% x1    --  sum_k (E(x_1))
% Q1    --  sum_k (E(x_1 * x_1'))
% K     --  number of trials
%
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
function [x0, Q0] = M_step_init (x1, Q1, K)

        x0     = x1/K;
        Q0     = Q1/K - x0*x0';