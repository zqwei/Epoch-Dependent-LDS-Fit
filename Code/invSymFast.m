% 
% inv(A) where A is positive-definite and symmetric.
% 
% @ 2014 Ziqiang Wei
% weiz@janelia.hhmi.org
%
function invA = invSymFast(A)
    if ~isequal(A, A'); A = (A+A')/2; end;
%     [invA, is_success] = invSym(A);
%     if ~is_success && nargout == 1; invA = inv(A); end 
    
    invA = A\eye(size(A));