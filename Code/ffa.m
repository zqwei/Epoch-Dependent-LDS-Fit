%
% [L,Ph,LL]=ffa(X,K,cyc,tol);
% 
% Fast Maximum Likelihood Factor Analysis using EM
%
% X - data matrix: N x D: N: trials; D: xDim
% K - number of factors
% cyc - maximum number of cycles of EM (default 100)
% tol - termination tolerance (prop change in likelihood) (default 0.0001)
%
% L - factor loadings 
% Ph - diagonal uniquenesses matrix
% LL - log likelihood curve
%
% Iterates until a proportional change < tol in the log likelihood 
% or cyc steps of EM 
%
% @ 1996 Zoubin Ghahramani
% 
% Modified by
% @ 2014 Ziqiang Wei
% weiz@janelia.hhmi.org
%
% Adding invSymFast to increase the speed of computing inv_symmetric
% matrices;
%
% Other optimization of previous code.
%

function [L, Ph] = ffa(X, K, cyc, tol)

    if nargin<4; tol = 0.0001; end;
    if nargin<3; cyc = 100; end;

    [N, D] = size(X);

    mean_X = mean(X);
    X      = X - ones(N,1) * mean_X;
    XX     = X'*X/N;
    diagXX = diag(XX);

    rng('shuffle');
    cX     = cov(X);
    scale  = det(cX)^(1/D);
    if scale == 0
        scale = 1e-8;
    end
    % initializing L with a randn matrix
    L      = randn(D,K)*sqrt(scale/K);
    Ph     = diag(cX);
    I      = eye(K);
    lik    = 0;
    const  = -D/2*log(2*pi);

    for n_cyc = 1:cyc
      
      %%%%%%%%%%%%%%%%
      %%%% E Step %%%%
      %%%%%%%%%%%%%%%%
      
      oldlik = lik;
      MM     = invSymFast(diag(Ph) + L*L');
%       if ~is_success; 
%           return; 
%       end
      
      beta   = L'*MM;
      XXbeta = XX*beta';
%       EZ     = X * beta';
      EZZ    = I-beta*L +beta*XXbeta;

      %%%% ---- log likelihood ----- %%%%
      lik  = N*const + 0.5*N*logdet(MM) - 0.5*N*trace(MM*XX);
      
      if n_cyc <= 2    
        likbase = lik;
      elseif lik < oldlik     
        disp('VIOLATION');
      elseif (lik-likbase) < (1+tol)*(oldlik-likbase) ||~isfinite(lik) 
        fprintf('FFA finished at cycle %i lik %g \n',n_cyc,lik);
        break;
      end;

      %%%%%%%%%%%%%%%%
      %%%% M Step %%%%
      %%%%%%%%%%%%%%%%
      L    = XXbeta * invSymFast(EZZ);
      Ph   = diagXX - diag(L * XXbeta');
      
    end % for