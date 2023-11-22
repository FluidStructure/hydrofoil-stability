function [y]=ferbsim(M,K,F,g1,g2,g3,EI,h,l0,thf,tf)
%----------------------------------------------------------------------------
%  Purpose:
%     This MATLAB m-file ferbsim.m is a simulation program for
%     a rotating flexible beam attached to a base. The mathematicla model
%     is created from ferobem.m as system mass and stifness matrices
%
%  Synopsis:
%     [y] = frobsim(M,K,F,g1,g2,g3,EI,h,l0,thf,tf)  
%
%  Variable Description:
%     Input parameters - M, K, F - System matrices
%                        g1, g2, g3 - Feedback gains	
%                        EI, h, l0 - Parameters for boundary force calculation
%                        thf, tf - Final angle and final simulation time
%     Output parameter - y 
%-----------------------------------------------------------------------------

F=zeros(9,1);
F(1,1)=1;EI=1.1118*10^4;h=47.57/4;l0=3.5;thf=1;g1=100;g2=200;g3=0.0;tf=40;

[n,n]=size(M);
I=eye(n);

%  Build closed-loop system matrices   

K(1,1)=g1;
K1I=EI*(l0*[-12/h^3,6/h^2]-[6/h^2,-2/h]);

K(1,2:3)=K(1,2:3)+g3*K1I;


Damp=0*I;
Damp(1,1)=g2;

%  Generate state space form for simulation purpose

A=[0*I,I;-inv(M)*K,-inv(M)*Damp];
B=[zeros(n,1);inv(M)*F];
C=eye(2*n);
D=zeros(2*n,1);

%  Perform simulation using felresp.m function 

nstep=500;
dt=tf/nstep;

t=0:dt:tf-dt;

u= g1*thf*ones(500,1);

x0=zeros(2*n,1);
y=felresp(A,B,C,D,u,t,x0);

%----------------------------------------------------------------------