function [dhdx,dhdy,dhdz]=federiv3(nnel,dhdr,dhds,dhdt,invjacob)

%------------------------------------------------------------------------
%  Purpose:
%     determine derivatives of 3-D isoparametric shape functions with 
%     respect to physical coordinate system
%
%  Synopsis:
%     [dhdx,dhdy,dhdz]=federiv3(nnel,dhdr,dhds,dhdt,invjacob)      
%
%  Variable Description:
%     dhdx - derivative of shape function w.r.t. physical coordinate x
%     dhdy - derivative of shape function w.r.t. physical coordinate y
%     dhdz - derivative of shape function w.r.t. physical coordinate z
%     nnel - number of nodes per element   
%     dhdr - derivative of shape functions w.r.t. natural coordinate r
%     dhds - derivative of shape functions w.r.t. natural coordinate s
%     dhdt - derivative of shape functions w.r.t. natural coordinate t
%     invjacob - inverse of 3-D Jacobian matrix
%------------------------------------------------------------------------

 for i=1:nnel
 dhdx(i)=invjacob(1,1)*dhdr(i)+invjacob(1,2)*dhds(i)+invjacob(1,3)*dhdt(i);
 dhdy(i)=invjacob(2,1)*dhdr(i)+invjacob(2,2)*dhds(i)+invjacob(2,3)*dhdt(i);
 dhdz(i)=invjacob(3,1)*dhdr(i)+invjacob(3,2)*dhds(i)+invjacob(3,3)*dhdt(i);
 end
