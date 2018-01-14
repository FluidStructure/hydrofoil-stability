%----------------------------------------------------------------------
%
%  This program {\it compen.m} demonstrates a dynamic compensator 
%  design and diplay the simulation result. A finte beam element is
%  adopted as a system. For obeserver gain design, the LQR technique
%  is used, and the dynamic siumlation is performed using a MATLAB
%  control toolbox routine {\it lsim.m}.
%
%---------------------------------------------------------------------

clear
% Provide the system mass and stiffness matrix

M =[0.5571         0    0.0964   -0.3482
         0    3.2143    0.3482   -1.2054
    0.0964    0.3482    0.2786   -0.5893
   -0.3482   -1.2054   -0.5893    1.6071];

K =[2.4178         0   -1.2089    9.0667
         0  181.3333   -9.0667   45.3333
   -1.2089   -9.0667    1.2089   -9.0667
    9.0667   45.3333   -9.0667   90.6667];

F=[1;0;0;0];

% Transform into the first order state space form equation

A=[0*eye(4),eye(4);-inv(M)*K,0*eye(4)];

B=[0*ones(4,1);inv(M)*F];

C=[1 0 0 0 0 0 0 0];

%  Use the {\it flqr.m} function to design the observer gain and the
%  full state feedback gain

[G,Sc]=felqr(A, B,  1000*eye(8), 0.01); 
[L,So]=felqr(A',C', eye(8),0.01);

% Now build the total closed-loop system for both closed-loop system and observer

Atot=[A, -B*G;L'*C, A-L'*C-B*G];
Btot=[B;B];
Ctot=eye(16);
Dtot=eye(16,1);

%  Define simulation time, control input, and initial conditions

t=0:0.03:6.0-0.03;
u=zeros(200,1);
x0=zeros(16,1);
x0(1,1)=0.1; x0(2,1)=-0.3; x0(3,1)=0.2;

%  Use MATLAB function {\it felresp.m} to simulate the total system

[x,y]=felresp(Atot,Btot,Ctot,Dtot,x0,u,t);

%  Plot the result

subplot(221)
plot(t,y(:,1),'-',t,y(:,9),':');
xlabel('Time(sec)');
ylabel('Displacement');

subplot(222)
plot(t,y(:,2),'-',t,y(:,10),':');
xlabel('Time(sec)');
ylabel('Rotation');

subplot(223)
plot(t,y(:,5),'-',t,y(:,13),':');
xlabel('Time(sec)');
ylabel('Linear velocity');

subplot(224)
plot(t,y(:,6),'-',t,y(:,14),':');
xlabel('Time(sec)');
ylabel('Angular velocity');

