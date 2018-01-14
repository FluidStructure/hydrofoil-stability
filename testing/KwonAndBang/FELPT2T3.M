function [m]=felpt2t3(x1,y1,x2,y2,x3,y3)

%-------------------------------------------------------------------
%  Purpose:
%     element matrix for transient term of two-dimensional 
%     Laplace's equation using linear triangular element
%
%  Synopsis:
%     [m]=felpt2t3(x1,y1,x2,y2,x3,y3) 
%
%  Variable Description:
%     m - element stiffness matrix (size of 3x3)   
%     x1, y1 - x and y coordinate values of the first node of element
%     x2, y2 - x and y coordinate values of the second node of element
%     x3, y3 - x and y coordinate values of the third node of element
%-------------------------------------------------------------------

% element matrix

 A=0.5*(x2*y3+x1*y2+x3*y1-x2*y1-x1*y3-x3*y2); % area of the triangule
 
 m = (A/12)* [ 2  1   1;
               1  2   1;
               1  1   2 ];
   

