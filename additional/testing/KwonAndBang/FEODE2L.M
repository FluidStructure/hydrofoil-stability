function [k]=feode2l(acoef,bcoef,ccoef,eleng)

%-------------------------------------------------------------------
%  Purpose:
%     element matrix for (a u'' + b u' + c u)
%     using linear element
%
%  Synopsis:
%     [k]=feode2l(acoef,bcoef,ccoef,eleng) 
%
%  Variable Description:
%     k - element matrix (size of 2x2)   
%     acoef - coefficient of the second order derivative term 
%     bcoef - coefficient of the first order derivative term
%     ccoef - coefficient of the zero-th order derivative term
%     eleng - element length
%-------------------------------------------------------------------

% element matrix

 a1=-(acoef/eleng);  a2=bcoef/2;  a3=ccoef*eleng/6;
 k=[ a1-a2+2*a3   -a1+a2+a3;...
     -a1-a2+a3    a1+a2+2*a3];
   

