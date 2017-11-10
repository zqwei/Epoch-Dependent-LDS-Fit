function   [X,Y] = SimulateLDS_MultiStage(LDS, T, nTrial, timePoints)

% Simulate a Linear Dynamical System
%
%
% INPUT:
%  
%  LDS is a Linear Dynamical System model with fields:
%
%    A,C,Q,R,x0,V0
%
%  Corresponding to the Model: 
%
%    x(t+1) = A*x(t) + q(t)
%    y(t)   = C*x(t) + r(t)
%
%    cov([q,r])=[Q 0; 0 R]
%
%    x0    initial state estimate as column vector - (nx)x1 
%          (or a cell array with a vector for each "run")
%    V0    cov of the initial state estimate - (nx)x(nx)
%
% 
%  T     Number of time steps
%        for one experiment:  
%           T is scalar
%        for one or more experiments:
%           T is vector of length E containing scalar trials counts
%
%
% OUTPUT
%
%  X    State
%  Y    Output
%
%  For one experiment:  
%     X - (nx)xT
%     Y - (ny)xT
%  For more than one experiment:
%     X,Y are cell arrays of length E containing
%       X - (nx)xT(e)
%       Y - (ny)xT(e)

% Copyright (C) 2005 Philip N. Sabes
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
%
%
% Modified by Ziqiang Wei
% weiz@janelia.hhmi.org 
%
% 07/22/2014

%%%%%%%%%%%  Dimensions
ny   = size(LDS.C,1);
nx   = size(LDS.A,1);   % k= dim of states

X    = zeros(nx, T, nTrial);
Y    = zeros(ny, T, nTrial);
sqV  = sqrtm(LDS.V0);
x    = LDS.x0(:,ones(nTrial,1)) + sqV*randn(nx,nTrial);

for t = 1:T
    nStage   = sum(t>timePoints)+1;
    sqQ      = sqrtm(LDS.Q(:,:,nStage));
    sqR      = sqrtm(LDS.R(:,:,nStage));
    A        = LDS.A(:,:,nStage);
    C        = LDS.C(:,:,nStage);
    y        = C*x + sqR*randn(ny,nTrial);  % output eq
    X(:,t,:) = x;
    Y(:,t,:) = y;
    x        = A*x + sqQ*randn(nx,nTrial);
end
