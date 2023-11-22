function [kinmtps]=fekineps(nnel,dhdx,dhdy,shape)

%------------------------------------------------------------------------
%  Purpose:
%     determine the kinematic matrix expression relating shear strains 
%     to rotations and displacements for shear deformable plate bending
%
%  Synopsis:
%     [kinmtps]=fekineps(nnel,dhdx,dhdy,shape) 
%
%  Variable Description:
%     nnel - number of nodes per element
%     dhdx - derivatives of shape functions with respect to x   
%     dhdy - derivatives of shape functions with respect to y
%     shape - shape function
%------------------------------------------------------------------------

 for i=1:nnel
 i1=(i-1)*3+1;  
 i2=i1+1;
 i3=i2+1;
 kinmtps(1,i1)=-shape(i);
 kinmtps(1,i3)=dhdx(i);
 kinmtps(2,i2)=-shape(i);
 kinmtps(2,i3)=dhdy(i);
 end
