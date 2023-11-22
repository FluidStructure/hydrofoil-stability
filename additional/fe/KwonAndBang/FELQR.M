function [G,S]=felqr(A,B,Q,R);
%------------------------------------------------------------------------------------------
%  Purpose:
%     The function subroutine felqr.m calculates the feedback gain matrix by
%     Linear Quadratic Regulator(LQR) technique. The given system is
%  
%                         xdot = Ax + Bu,     u = -Gx
%
%     and the cost function to be minimized is defined as
%
%                    J=(1/2)integral(x^TQx+u^TRu)dt
%  Synopsis:
%     [G,S]=flqr(A,B,Q,R)
%
%  Variable Description:
%     Input arguments   - A, B, Q, R
%     Output parameters - G = R^{-1}G^TS : feedback gain matrix
%                        S : Solution of the Algebraic Ricatti Equation (ARE)
%                                       AS+A^TS-SBR^{-1}S+Q=0
%  Notes:
%     i) (A,B) should be controllable.
%     ii) Q is at least positive semi-definite.
%         R is at least positive definite.
%------------------------------------------------------------------------------------------


H=[A -B*inv(R)*B';                 	             %  Build the Hamiltonian matrix
  -Q        -A'];

[V,D]=eig(H);					             %  Solve eigenvalue problem

n=size(A); twon=max(size(H));


% Normalized each eigenvector to unity magnitude

av=abs(V);
magav=av'*av;

dmagav=diag(magav);

V=V*sqrt(inv(diag(dmagav)));					%  Normalize the eigenvector
%-------------------------------------------------------------------------------------------
%  Sort the eigenvalues with stable real parts
%-------------------------------------------------------------------------------------------

rel=real(diag(D));

nindex=[];pindex=[];

for i=1:twon

if(rel(i)<=0)
nindex=[nindex i];
else
pindex=[pindex i];
end

end

V=V(:, [pindex,nindex]);					% Rearrange the eigenvector order

%------------------------------------------------------------------------------------------
%    Compute the feedback gain matrix and Riccati Matrix
%------------------------------------------------------------------------------------------

S=real(V(n+1:twon,n+1:twon)*inv(V(1:n,n+1:twon)));

G=real(inv(R)*B'*S);
%------------------------------------------------------------------------------------------ 
