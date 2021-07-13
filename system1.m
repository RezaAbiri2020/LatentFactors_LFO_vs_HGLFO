function [ xdot ] = system1( t,x)

load('AMatrix.mat')

if nargin==1
    x=t;
end
xdot(1,1)=A(1,1)*x(1)+A(1,2)*x(2);
xdot(2,1)=A(2,1)*x(1)+A(2,2)*x(2);

end