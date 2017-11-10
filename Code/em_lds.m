%
% Ph = em_lds(Yk_t, m, cyc, tol, C, R, A, Q, x0, Q0)
%
% Parameter estimation for linear dynamic system with external inputs
%
% Inputs:
%
% Y   -- n x (T x K) observation data matrix
% K   -- number of trials
% n   -- dimension of observation
% T   -- length of obervations in each trial 
% m   -- dimension of state variable
% cyc -- maximum number of cycles of EM (default: 1000)
% tol -- termination tolerance (% change in likelihood) (default: 0.01%)
%
% Output Ph (struct) ->
% 
% Ph.A  -- state transition matrix
% Ph.Q  -- state noise covarince
% Ph.x0 -- initial state mean
% Ph.Q0 -- initial state covariance
% Ph.C  -- observation matrix
% Ph.R  -- observation covariance
% Ph.LL -- log likelihood
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
% This is main function that handles EM iterations until change of LL < tol
% or cyc step reaches maximum cyc value.
%
% Ver: 2.1 adding individual m-step functions
% Ver: 2.0
% Ver: 1.0 renamed as lds_old
%
% @ 2014 Ziqiang Wei
% weiz@janelia.hhmi.org
% 
%




function Ph = em_lds(Yk_t, m, cyc, tol, C, R, A, Q, x0, Q0, timePoint)
    
    [yDim, T, K] = size(Yk_t);

    LL        = nan(cyc,1);
    nt        = size(timePoint,2) - 1;
        
    Q_yy_tt   = nan(yDim, yDim, nt);
    
    for nt_now = 1 : nt
        Y_all               = ...
            reshape(Yk_t(:,(timePoint(nt_now)+1):timePoint(nt_now + 1),:), yDim, []);
        
        Q_yy_tt(:,:,nt_now) = Y_all * Y_all';
        
    end
    
    Xk_t  = nan(m, T, K);
   
    disp ('----------------------------------------------');
    disp ('------- Starting EM Iteration for LDS Fit ----');
    disp ('----------------------------------------------');
    
    is_converage = false;
    
    for n_cyc = 1:cyc % EM iteration
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% ---- E step
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
        Xk_t_old   = Xk_t;
        
        Est        = kalman_forward_backward(Yk_t, A, Q, C, R, x0, Q0, timePoint);
        X_0        = Est.X_0;
        Q_yx_tt    = Est.Q_yx_tt;
        Q_xx_tt    = Est.Q_xx_tt;
        Q_xx_ts    = Est.Q_xx_ts;
        Q_xx_ss    = Est.Q_xx_ss;
        Q_xx_uu    = Est.Q_xx_uu;
        Xk_t       = Est.Xk_t;
        lik        = Est.lik;

        LL(n_cyc)  = lik;
        
        fprintf('cycle %g lik %g',n_cyc,LL(n_cyc)); 
        
        % If-Loop for "end of EM"
        if n_cyc <= 2
            likbase = lik; 
        elseif lik < oldlik
           fprintf('Violation'); 
           if n_cyc == cyc
               error('Violation and a new fit would start'); %break;
           end
        elseif (lik - likbase) < (1 + tol) * (oldlik-likbase) || ~isfinite(lik)
            if n_cyc < cyc; LL = LL(1:n_cyc); end
            is_converage = true;
            fprintf('\n'); 
            disp   ('------------ End of LDS EM -------------------');
            disp   ('----------------------------------------------');
            break;
        else
            Xk_t_old   = Xk_t;
        end        
        fprintf('\n');
        
        oldlik = lik;
        
%         figure;plot(squeeze(Xk_t))
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% ---- M step
        %%%      estimation of Ph
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        x0       = mean(X_0,2);
        Q0       = squeeze(Q_xx_tt(:,:,1) - Q_xx_uu(:,:,1))/K - x0*x0';
        
        for nt_now = 1: nt            
            length_T  = K*(timePoint(nt_now + 1) - timePoint(nt_now));            
            [C_now, R_now, ~] = M_step_obs (Q_yy_tt(:,:,nt_now), Q_yx_tt(:,:,nt_now), Q_xx_tt(:,:,nt_now), length_T);
            C(:,:,nt_now)     = C_now;
            R(:,:,nt_now)     = R_now;
            
            [A_now, Q_now, ~] = M_step_state (Q_xx_ts(:,:,nt_now), Q_xx_ss(:,:,nt_now), Q_xx_uu(:,:,nt_now), length_T - K);
            A(:,:,nt_now) = A_now;
            Q(:,:,nt_now) = Q_now;  
        end
        
    end
    
    
    % summarize outputs to structure Ph
    
    Ph.A  = A;
    Ph.Q  = Q;
    Ph.x0 = x0;
    Ph.Q0 = Q0;
    Ph.C  = C;
    Ph.R  = R;
    Ph.LL = LL;
    Ph.is_success = true;
    Ph.Xk_t = Xk_t_old;
    Ph.is_converage = is_converage;