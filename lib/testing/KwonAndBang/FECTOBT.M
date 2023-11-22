function [Ctobty,rrank,ccond]=fectobt(A,B)
%------------------------------------------------------------------------------------------
%  Purpose:
%     The function subroutine fectobt.m calculates controllability matrix and/or
%     oberservability of a system described in state space form
%
%                      xdot = Ax + Bu
%
%   i) For controllability test, the input argument should follow as
%
%                        fectobt(A,B)
%
%   ii) For observability test, we should provide the input argument as
%
%                        fectobt(A^T, C^T)  : ( )^T is transpose of ( )
%  Synopsis:
%     [Ctobty,rrank,ccond]=fectobt(A,B)
%
%  Variable Description:
%     Oputput parameters - Ctobty  : Controllability or observability matrix
%                          rrank  : rank of Ctobty which determine yes/no type answer
%                          ccond : Condition number of Ctobty
%-------------------------------------------------------------------------------------------


n=max(size(A));                 %  Find out the size of the system matrix

%------------------------------------------------------------------------------------------
%  Build the controllability/oberservability matrix (Ctobty)
%------------------------------------------------------------------------------------------

Ctobty=B;
Ao=A;

for i=1:n-1
Ctobty=[Ctobty Ao*B];
Ao=Ao*A;
end

rrank=rank(Ctobty);
ccond=cond(Ctobty);

%------------------------------------------------------------------------------------------
