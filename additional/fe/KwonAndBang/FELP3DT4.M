function [k]=felp3dt4(x,y,z)

%-------------------------------------------------------------------
%  Purpose:
%     element matrix for three-dimensional Laplace's equation
%     using four-node tetrahedral element
%
%  Synopsis:
%     [k]=felp3dt4(x,y,z) 
%
%  Variable Description:
%     k - element matrix (size of 4x4)   
%     x - x coordinate values of the four nodes
%     y - y coordinate values of the four nodes
%     z - z coordinate values of the four nodes
%-------------------------------------------------------------------

 xbar= [ 1  x(1)  y(1)  z(1);
         1  x(2)  y(2)  z(2);
         1  x(3)  y(3)  z(3);
         1  x(4)  y(4)  z(4) ];

 xinv = inv(xbar);

 vol = (1/6)*det(xbar);   % compute volume of tetrahedral

% element matrix

 k(1,1)=xinv(2,1)*xinv(2,1)+xinv(3,1)*xinv(3,1)+xinv(4,1)*xinv(4,1);
 k(1,2)=xinv(2,1)*xinv(2,2)+xinv(3,1)*xinv(3,2)+xinv(4,1)*xinv(4,2);
 k(1,3)=xinv(2,1)*xinv(2,3)+xinv(3,1)*xinv(3,3)+xinv(4,1)*xinv(4,3);
 k(1,4)=xinv(2,1)*xinv(2,4)+xinv(3,1)*xinv(3,4)+xinv(4,1)*xinv(4,4);
 k(2,1)=k(1,2);
 k(2,2)=xinv(2,2)*xinv(2,2)+xinv(3,2)*xinv(3,2)+xinv(4,2)*xinv(4,2);
 k(2,3)=xinv(2,2)*xinv(2,3)+xinv(3,2)*xinv(3,3)+xinv(4,2)*xinv(4,3);
 k(2,4)=xinv(2,2)*xinv(2,4)+xinv(3,2)*xinv(3,4)+xinv(4,2)*xinv(4,4);
 k(3,1)=k(1,3);
 k(3,2)=k(2,3);
 k(3,3)=xinv(2,3)*xinv(2,3)+xinv(3,3)*xinv(3,3)+xinv(4,3)*xinv(4,3);
 k(3,4)=xinv(2,3)*xinv(2,4)+xinv(3,3)*xinv(3,4)+xinv(4,3)*xinv(4,4);
 k(4,1)=k(1,4);
 k(4,2)=k(2,4);
 k(4,3)=k(3,4);
 k(4,4)=xinv(2,4)*xinv(2,4)+xinv(3,4)*xinv(3,4)+xinv(4,4)*xinv(4,4);
 k=vol*k;
