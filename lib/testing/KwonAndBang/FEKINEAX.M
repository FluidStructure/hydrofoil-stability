function [kinmtxax]=fekineax(nnel,dhdx,dhdy,shape,radist)

%------------------------------------------------------------------------
%  Purpose:
%     determine kinematic equations between strains and displacements
%     for axisymmetric solids
%
%  Synopsis:
%     [kinmtxax]=fekineax(nnel,dhdx,dhdy,shape,radist) 
%                        
%  Variable Description:
%     nnel - number of nodes per element
%     shape - shape functions
%     dhdx - derivatives of shape functions with respect to x   
%     dhdy - derivatives of shape functions with respect to y
%     radist - radial distance of integration point or central point
%              for hoop strain component
%------------------------------------------------------------------------

 for i=1:nnel
 i1=(i-1)*2+1;  
 i2=i1+1;
 kinmtxax(1,i1)=dhdx(i);
 kinmtxax(2,i1)=shape(i)/radist;
 kinmtxax(3,i2)=dhdy(i);
 kinmtxax(4,i1)=dhdy(i);
 kinmtxax(4,i2)=dhdx(i);
 end
