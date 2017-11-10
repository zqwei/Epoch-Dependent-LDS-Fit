% This code is a copy of Stochastic Realization
% Appendix D3
% Subspace Mehtods for System Idenfication
% T. Katayama, 2005

function [A, C, Cb, K, R] = stochastic(yk, n, k)
    
    [p, Ndat, totTrial] = size(yk);
    N          = Ndat - 2*k;
    
    Rfp        = 0;
    RR         = 0;
    for nTrial = 1:totTrial
        y      = yk(:,:,nTrial);
        ii        = 0;

        for i     = 1:p:2*k*p-p+1
            ii    = ii + 1;
            Y(i:i+p-1,:) = y(:,ii:ii+N-1);
        end

        Ypp       = Y(1:k*p,:);
        for i     = 1:k
            j     = (k-1)*p+1;
            Yp(j:j+p-1,:) = Ypp((i-1)*p+1:i*p,:); % Yp := Y check
        end

        Yf        = Y(k*p+1:2*k*p,:);
        Rfp       = Rfp + (Yf*Yp')/N;
        RR        = RR  + (Yf*Yf')/N;
    end
    
    Rfp           = Rfp / nTrial;
    RR            = RR  / nTrial;
    
    %%% Plot of SVD of Rfp to show the possible largest eigenvalues
    
    plot(svd(Rfp),'ok');
    
    [U, S, V]     = svd(Rfp);
    S2            = sqrtm(S(1:n,1:n));
    Ok            = U(:,1:n) * S2;
    Ck            = S2*V(:,1:n)';
    A             = Ok(1:k*p-p,:)\Ok(p+1:k*p,:);
    C             = Ok(1:p,:);
    Cb            = Ck(1:n,1:p)';
    R0            = RR(1:p,1:p);
    [P, ~, G, ~]  = dare(A', C', zeros(n,n), -R0, Cb');
    K             = G';
    R             = R0 - C*P*C';