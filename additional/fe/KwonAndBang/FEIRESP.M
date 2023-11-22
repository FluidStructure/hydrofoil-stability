function [eta,yim]=feiresp(M,K,F,u,t,C,q0,dq0);

%------------------------------------------------------------------------
%  Purpose:
%     The function subroutine feiresp.m calulates impulse response 
%     for a given structural system using modal analysis. It uses modal 
%     coordinate equations to evaluate modal responses anaytically, then
%     convert modal coordinates into physical responses
%
%  Synopsis:
%     [eta,yim]=feiresp(M,K,F,u,t,C,q0,dq0)
%
%  Variable Description:
%     Input parameters  - M, K - Mass and stiffness matrices
%                         F - Input or forcing influence matrix
%                         u - Index for excitation
%                         t - Time of evaluation
%                         C - Output matrix
%                         q0, dq0 - Initial conditions
%     Output parameters - eta - modal coordinate response
%                         yim - physical coordinate response
%------------------------------------------------------------------------
disp('  ')
disp('Please wait!! - The job is being performed.')
%---------------------------------------------------------------
%  Solve the eigenvalue problem and normalized the eigenvectors
% --------------------------------------------------------------

[n,n]=size(M);[n,m]=size(F);
nstep=size(t');

[V,D]=eig(K,M);

[lambda,k]=sort(diag(D));  % Sort the eigenvaules and eigenvectors in 
                           % ascending order       
V=V(:,k); 

Factor=diag(V'*M*V);

Vnorm=V*inv(sqrt(diag(Factor)));      %  Eigenvectors are normailzed

omega=diag(sqrt(Vnorm'*K*Vnorm));  % Natural frequency in ascending order
 		   
Fnorm=Vnorm'*F;

%-------------------------------------------------------------------
%  Find out impulse response of each modal coordinate analytically
%-------------------------------------------------------------------

eta0=Vnorm'*M*q0;            % Initial conditions for modal coordinates
deta0=Vnorm'*M*dq0;          % - both displacement and velocity
eta=zeros(nstep,n);

for i=1:n               %  The responses are obtained for n modes
phase=omega(i)*t;
eta(:,i)=eta0(i)*cos(phase')+deta0(i)*sin(phase')/omega(i)+...
sin(phase')*Fnorm(i,u);
end

%-----------------------------------------------------------------------
%  Convert modal coordinate responses to physical coordinate responses
%-----------------------------------------------------------------------

yim=C*Vnorm*eta';

%-----------------------------------------------------------------------