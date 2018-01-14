function [m]=felpt2r4(xleng,yleng)

%-------------------------------------------------------------------
%  Purpose:
%     element matrix of transient term for two-dimensional Laplace's 
%     equation using four-node bilinear rectangular element
%
%  Synopsis:
%     [m]=felpt2r4(xleng,yleng) 
%
%  Variable Description:
%     m - element stiffness matrix (size of 4x4)   
%     xleng - element size in the x-axis
%     yleng - element size in the y-axis
%-------------------------------------------------------------------

% element matrix

  m=(xleng*yleng/36)*[4  2  1  2;     
                      2  4  2  1;
                      1  2  4  2;
                      2  1  2  4];

