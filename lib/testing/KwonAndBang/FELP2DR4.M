function [k]=felp2dr4(xleng,yleng)

%-------------------------------------------------------------------
%  Purpose:
%     element matrix for two-dimensional Laplace's equation
%     using four-node bilinear rectangular element
%
%  Synopsis:
%     [k]=felp2dr4(xleng,yleng) 
%
%  Variable Description:
%     k - element stiffness matrix (size of 4x4)   
%     xleng - element size in the x-axis
%     yleng - element size in the y-axis
%-------------------------------------------------------------------

% element matrix

 k(1,1)=(xleng*xleng+yleng*yleng)/(3*xleng*yleng);
 k(1,2)=(xleng*xleng-2*yleng*yleng)/(6*xleng*yleng);
 k(1,3)=-0.5*k(1,1);
 k(1,4)=(yleng*yleng-2*xleng*xleng)/(6*xleng*yleng);
 k(2,1)=k(1,2);   k(2,2)=k(1,1);   k(2,3)=k(1,4);   k(2,4)=k(1,3);
 k(3,1)=k(1,3);   k(3,2)=k(2,3);   k(3,3)=k(1,1);   k(3,4)=k(1,2);
 k(4,1)=k(1,4);   k(4,2)=k(2,4);   k(4,3)=k(3,4);   k(4,4)=k(1,1);
   

