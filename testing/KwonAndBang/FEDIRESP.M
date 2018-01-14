function [eta,yim]=fediresp(M,K,F,u,t,C,q0,dq0,a,b);

%----------------------------------------------------------------------------------
%  Purpose:
%     The function subroutine fediresp.m calulates impulse response 
%     for a damped structural system using modal analysis. It uses modal 
%     coordinate equations to evaluate modal responses anaytically, then
%     convert modal coordinates into physical responses
%
%  Synopsis:
%     [eta,yim]=fediresp(M,K,F,u,t,C,q0,dq0,a,b)
%
%  Variable Description:
%     Input parameters - M, K - Mass and stiffness matrices
%                        F - Input or forcing influence matrix
%                        t - Time of evaluation
%                        u - Index for excitation
%                        C - Output matrix
%                        q0, dq0 - Initial conditions
%                        a, b - Parameters for proportional damping [C]=a[M]+b[K]
%      Outpur parameters - eta - modal coordinate response
%                          y - physical coordinate response
%---------------------------------------------------------------------------------

disp('  ')
disp('Please wait!! - The job is being performed.')
%---------------------------------------------------------------
%  Solve the eigenvalue problem and normalize the eigenvectors
% --------------------------------------------------------------

[n,n]=size(M);[n,m]=size(F);
nstep=size(t');

[V,D]=eig(K,M);

[lambda,k]=sort(diag(D));  % Sort the eigenvaules and eigenvectors 
                                 
V=V(:,k); 

Factor=diag(V'*M*V);

Vnorm=V*inv(sqrt(diag(Factor)));      %  Eigenvectors are normailzed

omega=diag(sqrt(Vnorm'*K*Vnorm));  % Natural frequencies
 		   
Fnorm=Vnorm'*F;

%----------------------------------------------------------------------------
%  Compute modal damping matrix from the proportional damping matrix
%----------------------------------------------------------------------------
Modamp=Vnorm'*(a*M+b*K)*Vnorm;
zeta=diag((1/2)*Modamp*inv(diag(omega)))

if (max(zeta) >= 1),
disp('Warning - Your maximu damping ratio is grater than or equal to 1')
disp('You have to reselect a and b ')
end
pause
disp('If you want to continue, type return key')
end
%----------------------------------------------------------------------
%  Find out impulse response of each modal coordinate analytically
%-------------------------------------------------------------------

eta0=Vnorm'*M*q0;            % Initial conditions for modal coordinates
deta0=Vnorm'*M*dq0;          % - both displacement and velocity
eta=zeros(nstep,n);

for i=1:n               %  Responses are obtained for n modes
omegad=omega(i)*sqrt(1-zeta(i)^2);
phase=omegad*t;
tcons=zeta(i)*omega(i)*t;
eta(:,i)=exp(-tcons)'.*(eta0(i)*(cos(phase')+zeta(i)/sqrt(1-zeta(i)^2)*...
sin(phase'))+deta0(i)*sin(phase')/omegad+sin(phase')*Fnorm(i,u)/omegad);
end

%-----------------------------------------------------------------------
%  Convert modal coordinate responses to physical coordinate responses
%-----------------------------------------------------------------------

yim=C*Vnorm*eta';
%-----------------------------------------------------------------------