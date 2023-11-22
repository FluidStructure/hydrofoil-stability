function [kinmtss]=fekiness(nnel,dhdx,dhdy,shape)

%------------------------------------------------------------------------
%  Purpose:
%     determine the kinematic matrix expression relating shear strains 
%     to rotations and displacements for shear deformable plate bending
%
%  Synopsis:
%     [kinmtss]=fekiness(nnel,dhdx,dhdy,shape) 
%
%  Variable Description:
%     nnel - number of nodes per element
%     dhdx - derivatives of shape functions with respect to x   
%     dhdy - derivatives of shape functions with respect to y
%     shape - shape function
%------------------------------------------------------------------------

 for i=1:nnel
 i1=(i-1)*6+1;  
 i2=i1+1;
 i3=i2+1;
 i4=i3+1;
 i5=i4+1;
 i6=i5+1;
 kinmtss(1,i3)=dhdx(i);
 kinmtss(1,i5)=shape(i);
 kinmtss(2,i3)=dhdy(i);
 kinmtss(2,i4)=-shape(i);
 kinmtss(2,i6)=0.0;
 end
