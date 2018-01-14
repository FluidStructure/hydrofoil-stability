function [tr3d,xprime,yprime]=fetransh(xcoord,ycoord,zcoord,n)

%--------------------------------------------------------------
%  Purpose:
%     Compute direction cosines between three-dimensional
%    local and global coordinate axes
%
%  Synopsis:
%     [tr3d,xprime,yprime]=fetransh(xcoord,ycoord,zcoord,n) 
%
%  Variable Description:
%     xcoord - nodal x coordinates (4x1)   
%     ycoord - nodal y coordinates (4X1)
%     zcoord - nodal z coordinates (4X1)
%     n - number of nodes per element
%     tr3d - 3d transformation matrix from local to global axes
%     xprime - coordinate in terms of the local axes (4x1)
%     yprime - coordinate in terms of the global axes (4X1)
%
%  Note:
%     The local x-axis is defined in the direction from the first node
%     to the second node. Nodes 1, 2 and 4 define the local xy-plane.
%     The local z-axis is defined normal to the local xy-plane.
%     The local y-axis is defined normal to the x and z axes.
%--------------------------------------------------------------------------

%
%  compute direction cosines
%
v12x=xcoord(2)-xcoord(1); 
v12y=ycoord(2)-ycoord(1);
v12z=zcoord(2)-zcoord(1);
l12=sqrt(v12x^2+v12y^2+v12z^2);
v23x=xcoord(3)-xcoord(2);
v23y=ycoord(3)-ycoord(2);
v23z=zcoord(3)-zcoord(2);
l23=sqrt(v23x^2+v23y^2+v23z^2);
v34x=xcoord(4)-xcoord(3);
v34y=ycoord(4)-ycoord(3);
v34z=zcoord(4)-zcoord(3);
l34=sqrt(v34x^2+v34y^2+v34z^2);
v14x=xcoord(4)-xcoord(1);
v14y=ycoord(4)-ycoord(1);
v14z=zcoord(4)-zcoord(1);
l14=sqrt(v14x^2+v14y^2+v14z^2);
v13x=xcoord(3)-xcoord(1); 
v13y=ycoord(3)-ycoord(1);
v13z=zcoord(3)-zcoord(1);
l13=sqrt(v13x^2+v13y^2+v13z^2);
v1tx=v12y*v14z-v12z*v14y;
v1ty=v12z*v14x-v12x*v14z;
v1tz=v12x*v14y-v12y*v14x;
v1yx=v1ty*v12z-v1tz*v12y;
v1yy=v1tz*v12x-v1tx*v12z;
v1yz=v1tx*v12y-v1ty*v12x;
vxx=v12x/l12;
vxy=v12y/l12;
vxz=v12z/l12;
vyx=v1yx/sqrt(v1yx^2+v1yy^2+v1yz^2);
vyy=v1yy/sqrt(v1yx^2+v1yy^2+v1yz^2);
vyz=v1yz/sqrt(v1yx^2+v1yy^2+v1yz^2);
vzx=v1tx/sqrt(v1tx^2+v1ty^2+v1tz^2);
vzy=v1ty/sqrt(v1tx^2+v1ty^2+v1tz^2);
vzz=v1tz/sqrt(v1tx^2+v1ty^2+v1tz^2);

%
%   transformation matrix
% 
for i=1:2*n
i1=(i-1)*3+1;
i2=i1+1;
i3=i2+1;
tr3d(i1,i1)=vxx;
tr3d(i1,i2)=vxy;
tr3d(i1,i3)=vxz;
tr3d(i2,i1)=vyx;
tr3d(i2,i2)=vyy;
tr3d(i2,i3)=vyz;
tr3d(i3,i1)=vzx;
tr3d(i3,i2)=vzy;
tr3d(i3,i3)=vzz;
end

%  
%  compute nodal values in terms of local axes
%
alpa213=acos((l12^2+l13^2-l23^2)/(2*l12*l13));
alpa314=acos((l13^2+l14^2-l34^2)/(2*l13*l14));
alpa41y=2*atan(1)-alpa213-alpa314;
xprime(1)=0;  yprime(1)=0;
xprime(2)=l12;  yprime(2)=0;
xprime(3)=l13*cos(alpa213);  yprime(3)=l13*sin(alpa213);
xprime(4)=l14*sin(alpa41y);  yprime(4)=l14*cos(alpa41y);