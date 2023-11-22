function [x,y]=felresp(A,B,C,D,x0,u,t)
%---------------------------------------------------------------------------------
%  Purpose:
%     find the time response of a linear system driveny by initial condition 
%     and external input. The numerical algorithm used in this program is zero holder approximation
%     for control input for discretized system.
%
%  Synopsis:
%     [x,y]=felresp(A,B,C,D,x0,u,t)
%
%  Variable Description:
%     A, B, C, D; system matrices in
%
%     xdot = Ax + Bu,   y = Cx  + Du 
%
%     x0; initial condition vector for the state variables
%     t ; integration time at equal distance as t=0:dt:tf
%         dt- time step, tf - final time
%     u ; control input vector with as many rows as the size of t
%     x(y) ; state(output) vector
%
%  Notes:
%     The control input vector must have as many columns as the number of input
%---------------------------------------------------------------------------------

[n,m]=size(B);

%---------------------------------------------------------------
% Transform into discrete equation by zero-holder approximation
%---------------------------------------------------------------

Ts=t(2)-t(1);

Phi=expm(A*Ts);

Gamma=inv(A)*(Phi-eye(n))*B;

nc=max(size(t));

x=zeros(nc,n);
tx=zeros(n,1);

xi=x0;
tx=xi;

%   Calculate time responses x first

for i=1:nc
x(i,:)=tx';
tx=Phi*tx+Gamma*u(i,:)';
end

%  Calculate the output response by using y=Cx+Du

y=(C*x'+D*u')';
%-------------------------------------------------------------------------------------
