function [f]=fefx2l(xl,xr)

%-------------------------------------------------------------------
%  Purpose:
%     element vector for f(x)=x^2
%     using linear element
%
%  Synopsis:
%     [f]=fefx2l(xl,xr) 
%
%  Variable Description:
%     f - element vector (size of 2x1)   
%     xl - coordinate value of the left node 
%     xr - coordinate value of the right node
%-------------------------------------------------------------------

% element vector

 eleng=xr-xl;             % element length
 f=(1/(12*eleng))*[ (-4*xr*xl^3+xr^4+3*xl^4); ...
                    (-4*xr^3*xl+3*xr^4+xl^4)];
   

