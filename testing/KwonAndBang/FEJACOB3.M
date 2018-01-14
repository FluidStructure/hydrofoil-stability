function [jacob3]=fejacob3(nnel,dhdr,dhds,dhdt,xcoord,ycoord,zcoord)

%------------------------------------------------------------------------
%  Purpose:
%     determine the Jacobian for three-dimensional mapping
%
%  Synopsis:
%     [jacob3]=fejacob3(nnel,dhdr,dhds,dhdt,xcoord,ycoord,zcoord) 
%
%  Variable Description:
%     jacob3 - Jacobian for one-dimension
%     nnel - number of nodes per element   
%     dhdr - derivative of shape functions w.r.t. natural coordinate r
%     dhds - derivative of shape functions w.r.t. natural coordinate s
%     dhdt - derivative of shape functions w.r.t. natural coordinate t
%     xcoord - x axis coordinate values of nodes
%     ycoord - y axis coordinate values of nodes
%     zcoord - z axis coordinate values of nodes
%------------------------------------------------------------------------

 jacob3=zeros(3,3);

 for i=1:nnel
 jacob3(1,1)=jacob3(1,1)+dhdr(i)*xcoord(i);
 jacob3(1,2)=jacob3(1,2)+dhdr(i)*ycoord(i);
 jacob3(1,3)=jacob3(1,3)+dhdr(i)*zcoord(i);
 jacob3(2,1)=jacob3(2,1)+dhds(i)*xcoord(i);
 jacob3(2,2)=jacob3(2,2)+dhds(i)*ycoord(i);
 jacob3(2,3)=jacob3(2,3)+dhds(i)*zcoord(i);
 jacob3(3,1)=jacob3(3,1)+dhdt(i)*xcoord(i);
 jacob3(3,2)=jacob3(3,2)+dhdt(i)*ycoord(i);
 jacob3(3,3)=jacob3(3,3)+dhdt(i)*zcoord(i);
 end
