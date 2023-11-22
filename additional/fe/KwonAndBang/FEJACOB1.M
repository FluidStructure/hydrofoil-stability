function [jacob1]=fejacob1(nnel,dhdr,xcoord)

%-------------------------------------------------------------------
%  Purpose:
%     determine the Jacobian for one-dimensional mapping
%
%  Synopsis:
%     [jacob1]=fejacob1(nnel,dhdr,xcoord) 
%
%  Variable Description:
%     jacob1 - Jacobian for one-dimension
%     nnel - number of nodes per element   
%     dhdr - derivative of shape functions w.r.t. natural coordinate 
%     xcoord - x axis coordinate values of nodes
%-------------------------------------------------------------------

 jacob1=0.0;

 for i=1:nnel
 jacob1=jacob1+dhdr(i)*xcoord(i);
 end
