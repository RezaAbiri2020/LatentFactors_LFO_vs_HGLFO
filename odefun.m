function dxdt = odefun( t,x)

load('AMatrix.mat')
dxdt = zeros(2,1);
dxdt(1)=A(1,1)*x(1)+A(1,2)*x(2);
dxdt(2)=A(2,1)*x(1)+A(2,2)*x(2);

end