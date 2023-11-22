function [g]=febasgr(A,B,dc)
%------------------------------------------------------------------------------------------
%  Purpose:
%     The function subroutine febasgr.m caculates a feedbback gain for a single
%     input system by Bass-Gura formular.
% 
%  Synopsis:
%     [g]=febasgur(A,B,dc)
%
%     System equation :   xdot  = Ax + bu
%     
%  Variabe Description:
%     Input variables : dc - A vector consisting of desired closed-loop poles                    
%     Output :  g - A feedback gain vector.
%------------------------------------------------------------------------------------------


ao= poly(A);					% Calculate coefficient of the given system

alpha = poly(dc);                         % Calculate coefficient of the desired polynomial

[P,rank,cond]=fctobty(A,B);                     % Compute controllability matrix

n=max(size(A));

%------------------------------------------------------------------------------------------
%                 Build a Toeplitz matrix 
%------------------------------------------------------------------------------------------

Toep=zeros(n,n);

for i=1:n
Toep(i:n,i)=[ao(1:n-i+1)]';
end

g=[alpha(2:n+1)-ao(2:n+1)]*(inv(Toep))'*inv(P);                % Calculate the feedback gain 

g=real(g);					        % Take the real part of feedback gain 
%--------------------------------------------------------------------------------------------