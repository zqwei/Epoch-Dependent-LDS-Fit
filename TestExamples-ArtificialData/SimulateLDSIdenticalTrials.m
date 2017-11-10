function   [X, Y, Y_est] = SimulateLDSIdenticalTrials(LDS, T, nTrial)

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

%%%%%%%%%%% distribute LDS elements
FN = fieldnames(LDS);
for i=1:length(FN), 
  eval(sprintf('%s = LDS.%s;',FN{i},FN{i}));
end

%%%%%%%%%%%  Dimensions
E     = nTrial;
ny    = size(C,1);
nx    = size(A,1);   % k= dim of states

sqQ   = sqrtm(Q);
sqR   = sqrtm(R);
sqV   = sqrtm(V0);

X     = zeros(nx, T, nTrial);
Y     = zeros(ny, T, nTrial);
Y_est = zeros(ny, T, nTrial);

x     = LDS.x0(:,ones(nTrial,1)) + sqV*randn(nx,1)*ones(1,nTrial);
for t = 1:T
    y        = C*x + sqR*randn(ny,1)*ones(1,nTrial);  % output eq
    y_est    = C*x;
    X(:,t,:) = x;
    Y(:,t,:) = y;
    Y_est(:,t,:) = y_est;
    x        = A*x + sqQ*randn(nx,1)*ones(1,nTrial);
end
