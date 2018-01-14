function [kinmtpb]=fekinepb(nnel,dhdx,dhdy)

%--------------------------------------------------------------------------
%  Purpose:
%     determine the kinematic matrix expression relating bending curvatures 
%     to rotations and displacements for shear deformable plate bending
%
%  Synopsis:
%     [kinmtpb]=fekinepb(nnel,dhdx,dhdy) 
%
%  Variable Description:
%     nnel - number of nodes per element
%     dhdx - derivatives of shape functions with respect to x   
%     dhdy - derivatives of shape functions with respect to y
%--------------------------------------------------------------------------

 for i=1:nnel
 i1=(i-1)*3+1;  
 i2=i1+1;
 i3=i2+1;
 kinmtpb(1,i1)=dhdx(i);
 kinmtpb(2,i2)=dhdy(i);
 kinmtpb(3,i1)=dhdy(i);
 kinmtpb(3,i2)=dhdx(i);
 kinmtpb(3,i3)=0;
 end
