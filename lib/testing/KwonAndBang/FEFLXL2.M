function [k]=feflxl2(eleng)

%-------------------------------------------------------------------
%  Purpose:
%     element matrix for Cauchy-type boundary such as du/dn=a(u-b)
%     using linear element where a and b are known constants.
%
%  Synopsis:
%     [k]=feflxl2(eleng) 
%
%  Variable Description:
%     k - element vector (size of 2x2)   
%     eleng - length of element side with given flux
%-------------------------------------------------------------------

% element matrix

 k=(eleng/6)*[ 2  1;  
               1  2];
   

