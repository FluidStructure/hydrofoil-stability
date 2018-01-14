function [Omega,Phi,ModF]=femodal(M,K,F);
%-----------------------------------------------------------------------------------------------
%  Purpose:
%     The function subroutine femodal.m calulates modal parameters
%     for a given structural system. It calculates natural frequency and
%     eigenvector. The eigenvectors are normalized so that the modal
%     mass matrix becomes an indentity matrix.
%
%  Synopsis:
%     [Omega, Phi, ModF]=femodal(M,K,F)
%
%  Variable Description:
%     Input parameters - M, K - Mass and stiffness matrices
%                        F - Input or forcing influence matrix
%     Output parameters - Omega - Natural frequency(rad/sec) in ascending order
%                         Phi - Modal matrix with each column corresponding to the eigenvector.
%                         ModF - Modal input matrices.
%----------------------------------------------------------------------------------------------

disp('  ')
disp('Please wait!! - The job is being performed.')

%---------------------------------------------------------------
%  Solve the eigenvalue problem and normalized the eigenvectors
% --------------------------------------------------------------

[n,n]=size(M);[n,m]=size(F);

[V,D]=eig(K,M);

[lambda,k]=sort(diag(D));  % Sort the eigenvaules and eigenvectors in 
                           % ascending order       
V=V(:,k); 

Factor=diag(V'*M*V);

Phi=V*inv(sqrt(diag(Factor)));      %  Eigenvectors are normailzed

Omega=diag(sqrt(Phi'*K*Phi));  % Natural frequency in ascending order
 		   
ModF=Phi'*F;
%----------------------------------------------------------------------