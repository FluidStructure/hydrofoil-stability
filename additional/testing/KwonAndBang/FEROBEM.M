function [w,M,K]=ferobem(N,EI,rho,I_c,I_t,m_t,l_0,L)

%------------------------------------------------------------------------------------------
%  Purpose:
%     The MATLAB function subroutine ferobem.m produces a finite element modeling
%     of a rotating flexible beam attached to a rigid base.
%  
%  Synopsis: 
%      [w,M,K]=ferobem(N,EI,rho,I_c,I_t,m_t,l_0,L)
%
%  Variable Description:
%     Input  parameters: N    -  mumber of elements
%			 EI   -  elastic rigidity
%                        rho  -  linear mass density
%		         I_c(I_t) - moment of inertia of the center body(tip mass)
%			 m_t  - tip mass
%                        l_0 - radius of the center body
%                        L - beam length
%     Output: M, K - system mass, stiffness matrices
%	      w - natural frequency
%------------------------------------------------------------------------------------------

N=4;EI=1.1118*10^4;l_0=3.5;L=47.57;I_t=0.0018;m_t=0.156;I_c=9.06;rho=0.003;
%------------------------------------------------------------------------------------------
%          Calculate pure rigid body portion
%------------------------------------------------------------------------------------------

h=L/N;                                                          % Element length

M11t=I_t+m_t*(L+l_0)^2;M12t=[m_t*(l_0+L) I_t];M22t=[m_t 0; 0 I_t];M21t=M12t';

xi=0;
for i=1:N
M11r(1,i)=rho*h*((xi+l_0)^2+(xi+l_0+h)*(xi+l_0)+(xi+l_0+h)^2)/3;
xi=xi+h;
end

Mthth=M11t+sum(M11r)+I_c;

%------------------------------------------------------------------------------------------
%            Calculate element mass and stiffness matrices
%------------------------------------------------------------------------------------------

M33=rho*h*[156,-22*h;-22*h,4*h^2]/420;
M22=rho*h*[156, 22*h; 22*h,4*h^2]/420;
M23=rho*h*[54,-13*h;13*h,-3*h^2]/420;M32=M23';

K22=EI*[ 12, 6*h; 6*h,4*h^2]/h^3;
K23=EI*[-12, 6*h;-6*h,2*h^2]/h^3;K32=K23';
K33=EI*[ 12,-6*h;-6*h,4*h^2]/h^3;

%------------------------------------------------------------------------------------------
%   	     Calculate global mass and stiffness matrices
%------------------------------------------------------------------------------------------

Mqq(1:2,1:2)=M33+M22;Mqq(1:2,3:4)=M23;
for i=1:N-2
 NI=2*i+1;
 Mqq(NI:NI+1,NI-2:NI-1)=M32;Mqq(NI:NI+1,NI:NI+1)=M33+M22;
 Mqq(NI:NI+1,NI+2:NI+3)=M23;
end
Mqq(2*N-1:2*N,2*N-3:2*N-2)=M32; Mqq(2*N-1:2*N,2*N-1:2*N)=M33+M22t;


Kqq(1:2,1:2)=K33+K22;Kqq(1:2,3:4)=K23;
for i=1:N-2
 NI=2*i+1;
 Kqq(NI:NI+1,NI-2:NI-1)=K32;Kqq(NI:NI+1,NI:NI+1)=K33+K22;
 Kqq(NI:NI+1,NI+2:NI+3)=K23;
end
Kqq(2*N-1:2*N,2*N-3:2*N-2)=K32; Kqq(2*N-1:2*N,2*N-1:2*N)=K33;

%------------------------------------------------------------------------------------------
%  		Compute rigid and flexible coupling terms
%------------------------------------------------------------------------------------------

xi=0;
for i=1:N-1,    
M13=rho*h*[7*h/20+(xi+l_0)/2 -h^2/20-h*(xi+l_0)/12];
xi=xi+h;
M12=rho*h*[3*h/20+(xi+l_0)/2  h^2/30+h*(xi+l_0)/12];
Mthq(1,2*i-1:2*i)=M13+M12;
end
M13=rho*h*[7*h/20+(xi+l_0)/2 -h^2/20-h*(xi+l_0)/12];
Mthq(1,2*N-1:2*N)=M13+M12t;

%------------------------------------------------------------------------------------------
%  		Combine rigid, flexible, and coupling terms
%------------------------------------------------------------------------------------------

M=[Mthth,Mthq;Mthq',Mqq];                        % Global mass and stiffness matrices
K=[0,zeros(1,2*N);zeros(2*N,1),Kqq];

wo=sqrt(eig(K,M));                               % Compute natural frequencies
w=sort(wo);                    





      
