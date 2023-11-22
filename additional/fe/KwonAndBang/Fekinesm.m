function [kinmtsm]=fekinesm(nnel,dhdx,dhdy)

%------------------------------------------------------------------------
%  Purpose:
%     determine the kinematic equation between strains and displacements
%     for two-dimensional solids
%
%  Synopsis:
%     [kinmtsm]=fekinesm(nnel,dhdx,dhdy) 
%
%  Variable Description:
%     nnel - number of nodes per element
%     dhdx - derivatives of shape functions with respect to x   
%     dhdy - derivatives of shape functions with respect to y
%------------------------------------------------------------------------

 for i=1:nnel
 i1=(i-1)*6+1;  
 i2=i1+1;
 i3=i2+1;
 i4=i3+1;
 i5=i4+1;
 i6=i5+1;
 kinmtsm(1,i1)=dhdx(i);
 kinmtsm(2,i2)=dhdy(i);
 kinmtsm(3,i1)=dhdy(i);
 kinmtsm(3,i2)=dhdx(i);
 kinmtsm(3,i6)=0.0;
 end
