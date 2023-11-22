function [kinmtsb]=fekinesb(nnel,dhdx,dhdy)

%--------------------------------------------------------------------------
%  Purpose:
%     determine the kinematic matrix expression relating bending curvatures 
%     to rotations and displacements for shear deformable plate bending
%
%  Synopsis:
%     [kinmtsb]=fekinesb(nnel,dhdx,dhdy) 
%
%  Variable Description:
%     nnel - number of nodes per element
%     dhdx - derivatives of shape functions with respect to x   
%     dhdy - derivatives of shape functions with respect to y
%--------------------------------------------------------------------------

 for i=1:nnel
 i1=(i-1)*6+1;  
 i2=i1+1;
 i3=i2+1;
 i4=i3+1;
 i5=i4+1;
 i6=i5+1;
 kinmtsb(1,i5)=dhdx(i);
 kinmtsb(2,i4)=-dhdy(i);
 kinmtsb(3,i5)=dhdy(i);
 kinmtsb(3,i4)=-dhdx(i);
 kinmtsb(3,i6)=0;
 end
