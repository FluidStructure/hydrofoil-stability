function [shapet3,dhdrt3,dhdst3]=feisot3(rvalue,svalue)

%------------------------------------------------------------------------
%  Purpose:
%     compute isoparametric three-node triangular shape functions
%     and their derivatves at the selected (integration) point
%     in terms of the natural coordinate 
%
%  Synopsis:
%     [shapet3,dhdrt3,dhdst3]=feisot3(rvalue,svalue)  
%
%  Variable Description:
%     shapet3 - shape functions for three-node element
%     dhdrt3 - derivatives of the shape functions w.r.t. r
%     dhdst3 - derivatives of the shape functions w.r.t. s
%     rvalue - r coordinate value of the selected point   
%     svalue - s coordinate value of the selected point
%
%  Notes:
%     1st node at (0,0), 2nd node at (1,0), 3rd node at (0,1)
%------------------------------------------------------------------------

% shape functions

 shapet3(1)=1-rvalue-svalue;
 shapet3(2)=rvalue;
 shapet3(3)=svalue;

% derivatives

 dhdrt3(1)=-1;
 dhdrt3(2)=1;
 dhdrt3(3)=0;

 dhdst3(1)=-1;
 dhdst3(2)=0;
 dhdst3(3)=1;
