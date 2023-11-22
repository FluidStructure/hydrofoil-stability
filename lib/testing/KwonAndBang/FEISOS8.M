function [shapes8,dhdrs8,dhdss8,dhdts8]=feisos8(rvalue,svalue,tvalue)

%------------------------------------------------------------------------
%  Purpose:
%     compute isoparametric eight-node solid shape functions
%     and their derivatves at the selected (integration) point
%     in terms of the natural coordinate 
%
%  Synopsis:
%     [shapes8,dhdrs8,dhdss8,dhdts8]=feisos8(rvalue,svalue,tvalue)   
%
%  Variable Description:
%     shapes8 - shape functions for four-node element
%     dhdrs8 - derivatives of the shape functions w.r.t. r
%     dhdss8 - derivatives of the shape functions w.r.t. s
%     dhdts8 - derivatives of the shape functions w.r.t. t
%     rvalue - r coordinate value of the selected point   
%     svalue - s coordinate value of the selected point
%     tvalue - t coordinate value of the selected point
%
%  Notes:
%     1st node at (-1,-1,-1), 2nd node at (1,-1,-1)
%     3rd node at (1,1,-1), 4th node at (-1,1,-1)
%     5th node at (-1,-1,1), 6th node at (1,-1,1)
%     7th node at (1,1,1), 8th node at (-1,1,1)
%------------------------------------------------------------------------

% shape functions

 shapes8(1)=0.125*(1-rvalue)*(1-svalue)*(1-tvalue);
 shapes8(2)=0.125*(1+rvalue)*(1-svalue)*(1-tvalue);
 shapes8(3)=0.125*(1+rvalue)*(1+svalue)*(1-tvalue);
 shapes8(4)=0.125*(1-rvalue)*(1+svalue)*(1-tvalue);
 shapes8(5)=0.125*(1-rvalue)*(1-svalue)*(1+tvalue);
 shapes8(6)=0.125*(1+rvalue)*(1-svalue)*(1+tvalue);
 shapes8(7)=0.125*(1+rvalue)*(1+svalue)*(1+tvalue);
 shapes8(8)=0.125*(1-rvalue)*(1+svalue)*(1+tvalue);

% derivatives

 dhdrs8(1)=-0.125*(1-svalue)*(1-tvalue);
 dhdrs8(2)=0.125*(1-svalue)*(1-tvalue);
 dhdrs8(3)=0.125*(1+svalue)*(1-tvalue);
 dhdrs8(4)=-0.125*(1+svalue)*(1-tvalue);
 dhdrs8(5)=-0.125*(1-svalue)*(1+tvalue);
 dhdrs8(6)=0.125*(1-svalue)*(1+tvalue);
 dhdrs8(7)=0.125*(1+svalue)*(1+tvalue);
 dhdrs8(8)=-0.125*(1+svalue)*(1+tvalue);

 dhdss8(1)=-0.125*(1-rvalue)*(1-tvalue);
 dhdss8(2)=-0.125*(1+rvalue)*(1-tvalue);
 dhdss8(3)=0.125*(1+rvalue)*(1-tvalue);
 dhdss8(4)=0.125*(1-rvalue)*(1-tvalue);
 dhdss8(5)=-0.125*(1-rvalue)*(1+tvalue);
 dhdss8(6)=-0.125*(1+rvalue)*(1+tvalue);
 dhdss8(7)=0.125*(1+rvalue)*(1+tvalue);
 dhdss8(8)=0.125*(1-rvalue)*(1+tvalue);

 dhdts8(1)=-0.125*(1-rvalue)*(1-svalue);
 dhdts8(2)=-0.125*(1+rvalue)*(1-svalue);
 dhdts8(3)=-0.125*(1+rvalue)*(1+svalue);
 dhdts8(4)=-0.125*(1-rvalue)*(1+svalue);
 dhdts8(5)=0.125*(1-rvalue)*(1-svalue);
 dhdts8(6)=0.125*(1+rvalue)*(1-svalue);
 dhdts8(7)=0.125*(1+rvalue)*(1+svalue);
 dhdts8(8)=0.125*(1-rvalue)*(1+svalue);

