function S = smooth_kinematics(X,dim)
% function smooth_kinematics(X,dim)
% 
%   applies an sgolay smoother: Span=50, Degree=5
%   appies to each column of X, unless dim is specified
% 
%   Inputs:
%       X - matrix of pos/vel kinematics vs. time
%       *dim - dimension of X to act on (default=1)
% 
%   Outputs:
%       S - smoothed version of X

if ~exist('dim','var'), dim=1; end

if dim==2,
    X = X';
end

Ndim = size(X,2);
S = zeros(size(X));
for i=1:Ndim,
    S(:,i) = smooth(X(:,i),50,'sgolay',5);
end

if dim==2,
    S = S';
end
