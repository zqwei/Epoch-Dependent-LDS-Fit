% 
% function X = colvec(X)
% 
% transform a row vector to column vector
%
% @ 2014 Ziqiang Wei
% weiz@janelia.hhmi.org
%
%
function X = colvec(X)
    
    if isrow(X); X = X';end