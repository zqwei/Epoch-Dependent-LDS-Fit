% 
% log(det(A)) where A is positive-definite.
% This is faster and more stable than using log(det(A)).
% 
% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.
%
%
function y = logdet(A)


%     if ~all(eig(A)>0)
%         A = A + 1e-7*eye(size(A));
%     end
%     
%     A = (A + A')/2;
    
    U = chol(A);
    y = 2*sum(log(diag(U)));