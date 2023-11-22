function [m]=felpt3t4(x,y,z)

%-------------------------------------------------------------------
%  Purpose:
%     element matrix of transient term for three-dimensional Laplace's 
%     equation using four-node tetrahedral element
%
%  Synopsis:
%     [m]=felpt3t4(x,y,z) 
%
%  Variable Description:
%     m - element stiffness matrix (size of 4x4)   
%     xleng - element size in the x-axis
%     yleng - element size in the y-axis
%-------------------------------------------------------------------
  
  xbar = [ 1  x(1)  y(1)  z(1) ;
           1  x(2)  y(2)  z(2) ;
           1  x(3)  y(3)  z(3) ;
           1  x(4)  y(4)  z(4) ; ];
  vol = (1/6)*det(xbar);

% element matrix

  m=(vol/20)*[ 2  1  1  1;     
               1  2  1  1;
               1  1  2  1;
               1  1  1  2];

