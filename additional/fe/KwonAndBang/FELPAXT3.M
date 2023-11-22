function [k]=felpaxt3(r1,z1,r2,z2,r3,z3)

%-------------------------------------------------------------------
%  Purpose:
%     element matrix for axisymmetric Laplace equation
%     using three-node linear triangular element
%
%  Synopsis:
%     [k]=felpaxt3(r1,z1,r2,z2,r3,z3) 
%
%  Variable Description:
%     k - element stiffness matrix (size of 3x3)   
%     r1, z1 - r and z coordinate values of the first node of element
%     r2, z2 - r and z coordinate values of the second node of element
%     r3, z3 - r and z coordinate values of the third node of element
%-------------------------------------------------------------------

% element matrix

 A=0.5*(r2*z3+r1*z2+r3*z1-r2*z1-r1*z3-r3*z2); % area of the triangule
 rc=(r1+r2+r3)/3;    % r coordinate value of the centroid
 twopirc=8*atan(1)*rc;
 k(1,1)=((r3-r2)^2+(z2-z3)^2)/(4*A);
 k(1,2)=((r3-r2)*(r1-r3)+(z2-z3)*(z3-z1))/(4*A);
 k(1,3)=((r3-r2)*(r2-r1)+(z2-z3)*(z1-z2))/(4*A);
 k(2,1)=k(1,2);
 k(2,2)=((r1-r3)^2+(z3-z1)^2)/(4*A);
 k(2,3)=((r1-r3)*(r2-r1)+(z3-z1)*(z1-z2))/(4*A);
 k(3,1)=k(1,3);
 k(3,2)=k(2,3);
 k(3,3)=((r2-r1)^2+(z1-z2)^2)/(4*A);
 k=twopirc*k;  

