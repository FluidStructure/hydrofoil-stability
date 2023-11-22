function [k,m]=febeam2(el,xi,leng,sh,area,rho,ipt)

%--------------------------------------------------------------
%  Purpose:
%     Stiffness and mass matrices for C^0 beam element
%     nodal dof {v_1 theta_1 v_2 theta_2}
%
%  Synopsis:
%     [k,m]=febeam2(el,xi,leng,sh,area,rho,ipt) 
%
%  Variable Description:
%     k - element stiffness matrix (size of 4x4)     
%     m - element mass matrix (size of 4x4)
%     el - elastic modulus 
%     xi - second moment of inertia of cross-section
%     leng - length of the beam element
%     rho - mass density of the beam element (mass per unit volume)
%     sh - shear modulus
%     area - area of cross-section
%     ipt = 1: consistent mass matrix
%           2: lumped mass matrix
%           otherwise: diagonal mass matrix
%---------------------------------------------------------------

% stiffness matrix

 c=el*xi/leng;  
 d=(5/6)*sh*area/(4*leng);
 k= [ 4*d         2*d*leng     -4*d        2*d*leng;...
      2*d*leng    c+d*leng^2   -2*d*leng  -c+d*leng^2;...
     -4*d        -2*d*leng      4*d       -2*d*leng;...
      2*d*leng   -c+d*leng^2   -2*d*leng   c+d*leng^2];
   
% consistent mass matrix

 if ipt==1

    mm=rho*area*leng/6;
    m=mm*[2   0   1    0;...
          0   0   0    0;...
          1   0   2    0;...
          0   0   0    0];
  
% lumped or diagonal mass matrix

 else

    m=zeros(4,4);
    mass=rho*area*leng/2;
    m=mass*diag([1  0  1  0]);

end
 

