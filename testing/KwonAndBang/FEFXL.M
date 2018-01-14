function [f]=fefxl(xl,xr)

%-------------------------------------------------------------------
%  Purpose:
%     element vector for f(x)=x
%     using linear element
%
%  Synopsis:
%     [f]=fefxl(xl,xr) 
%
%  Variable Description:
%     f - element vector (size of 2x1)   
%     xl - coordinate value of the left node 
%     xr - coordinate value of the right node
%-------------------------------------------------------------------

% element vector

 eleng=xr-xl;             % element length
 f=[ eleng*(2*xl+xr)/6;  eleng*(xl+2*xr)/6];
   

