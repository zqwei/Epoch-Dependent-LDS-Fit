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

function Est = kalman_forward_backward(y, A_all, Q_all, C_all, R_all, x0, Q0, timePoint)

%%%%%%%% Initialized return result
    xDim          = size(A_all, 1);
    nt            = length(timePoint)-1;
    [yDim, T, K]  = size(y);   
    y             = permute(y, [1 3 2]); % yDim x K x T
    Est           = struct('X_0', zeros(xDim, K      ),...
                     'Xk_t',      zeros(xDim, T,    K),...
                     'Q_yx_tt',   zeros(yDim, xDim, nt),...
                     'Q_xx_tt',   zeros(xDim, xDim, nt),...
                     'Q_xx_ts',   zeros(xDim, xDim, nt),...
                     'Q_xx_ss',   zeros(xDim, xDim, nt),...
                     'Q_xx_uu',   zeros(xDim, xDim, nt),...
                     'is_success',false,...
                     'lik',0);

   
    %%%% Recording variables
   
    X_t         = zeros(xDim, K,    T);
    
    %%%% These variables are identical across trials
    Q_xx_tt     = zeros(xDim, xDim, T);
    Q_xx_ts     = zeros(xDim, xDim, T); 
    Q_yx_tt     = zeros(yDim, xDim, T);
    
    
    %%%% likelihood computed in each t step
    % See Notes for this equation
    sum_ll      = log(2*pi) * yDim * K * T;
    
    
    %%%% pre-computing of mutiplication matrices
    I_K         = eye(xDim);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%      Forward Pass
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q_xx_st        = zeros(xDim, xDim, T); 
    Q_xx_st(:,:,1) = Q0;
    Qst            = Q0;
    xst            = x0 * ones(1,K);

    for t = 1:T
        mt            = sum(t>timePoint);
        A             = A_all(:, :, mt);
        Q             = Q_all(:, :, mt);
        C             = C_all(:, :, mt);
        R             = R_all(:, :, mt);        
        yt            = squeeze(y(:, :,t)); % yDim * K
                
        
        ydiff         = yt - C*xst;
        MM            = R+C*Qst*C';
        invMM         = invSymFast(MM);     % yDim * yDim
        sum_ll        = sum_ll + K*logdet(MM) + sum(sum(ydiff.*(invMM*ydiff)));
        
        Kt            = Qst*C'*invMM;
        Qtt           = (I_K - Kt * C) * Qst;
        xtt           = xst + Kt*ydiff;

        % backward pass
        Q_xx_tt(:,:,t) = Qtt;
        X_t    (:,:,t) = xtt; % xDim x K x T

        if t<T
            % Next forward pass
            Qst               = A*Qtt*A' + Q;
            xst               = A*xtt;
            Q_xx_st(:,:,t+1)  = Qst;
        end
        
        

    end % Forward Pass

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%      Backward Pass
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Qu              = Qtt; % Q_T
    xu              = xtt; % X_T
    Qtu             = Qst; % Q_(T-1),T
    Qtt             = squeeze (Q_xx_tt(:,:,T-1));
    Qut             = A*Qtt - Kt * C*A*Qtt;
    Jt              = Qtt*A'*invSymFast(Qtu);
    yt              = squeeze(y(:, :,T));  % yDim * K
    Q_yx_tt(:,:,T)  = yt * xu';
    
    Q_xx_tt(:,:,T)  = xu * xu' + Qu*K;

    for t = (T-1):-1:1  
        
        mt            = sum(t>timePoint);
        A             = A_all(:, :, mt);
        
        xtt           = X_t(:,:,t);
        Qt            = Qtt + Jt*(Qu - Qtu)*Jt';
        xt            = xtt + Jt*(xu-A*xtt); 

        if t>1
            Qss     = squeeze (Q_xx_tt(:,:,t-1));
            Qst     = squeeze (Q_xx_st(:,:,t));                
            Js      = Qss*A'*invSymFast(Qst);
            Qts     = Qtt*Js' + Jt*(Qut - A*Qtt)*Js';
        end
                  
        Q_xx_tt(:,:,t)  = xt * xt' + Qt*K;
        X_t    (:,:,t)  = xt; % xDim x K x T
        yt              = squeeze(y(:, :,t));  % yDim * K
        Q_xx_ts(:,:,t)  = xu * xt' + Qut*K;
        Q_yx_tt(:,:,t)  = yt * xt';

        % Next backward pass
        Jt          = Js; 
        Qtt         = Qss;
        Qtu         = Qst;
        Qu          = Qt;
        Qut         = Qts;
        xu          = xt;

    end % Backward Pass
              
    Q_xx_ts(:,:,T) = zeros(xDim); 

    %%%%% Output
    
%     timePoint   = [0, timePoint, T];
    
    for nt_now  = 1:nt
        t_now   = (timePoint(nt_now) + 1): timePoint(nt_now + 1);
        Est.Q_xx_tt(:,:, nt_now)  = sum(Q_xx_tt(:,:,t_now),3);
        Est.Q_xx_ts(:,:, nt_now)  = sum(Q_xx_ts(:,:,t_now(1:end-1)),3);
        Est.Q_xx_ss(:,:, nt_now)  = sum(Q_xx_tt(:,:,t_now(1:end-1)),3);
        Est.Q_xx_uu(:,:, nt_now)  = sum(Q_xx_tt(:,:,t_now(2:end)),3);
        Est.Q_yx_tt(:,:, nt_now)  = sum(Q_yx_tt(:,:,t_now),3);
    end
    Est.X_0        = squeeze(X_t(:,:,1));
    Est.Xk_t       = permute(X_t, [1 3 2]);
    Est.is_success = true;
    Est.lik        = -sum_ll/2;