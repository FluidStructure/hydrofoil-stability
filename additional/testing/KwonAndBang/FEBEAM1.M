function [k,m]=febeam1(el,xi,leng,area,rho,ipt)

%-------------------------------------------------------------------
%  Purpose:
%     Stiffness and mass matrices for Hermitian beam element
%     nodal dof {v_1 theta_1 v_2 theta_2}
%
%  Synopsis:
%     [k,m]=febeam1(el,xi,leng,area,rho,ipt) 
%
%  Variable Description:
%     k - element stiffness matrix (size of 4x4)   
%     m - element mass matrix (size of 4x4)
%     el - elastic modulus 
%     xi - second moment of inertia of cross-section
%     leng - element length
%     area - area of beam cross-section
%     rho - mass density (mass per unit volume)
%     ipt = 1: consistent mass matrix
%           2: lumped mass matrix
%           otherwise: diagonal mass matrix
%-------------------------------------------------------------------

% stiffness matrix

 c=el*xi/(leng^3);
 k=c*[12      6*leng   -12       6*leng;...
      6*leng  4*leng^2 -6*leng   2*leng^2;...
     -12     -6*leng    12      -6*leng;...
      6*leng  2*leng^2 -6*leng   4*leng^2];
   
% consistent mass matrix

 if ipt==1

    mm=rho*area*leng/420;
    m=mm*[156      22*leng   54       -13*leng;...
          22*leng  4*leng^2  13*leng  -3*leng^2;...
          54       13*leng   156      -22*leng;...
         -13*leng -3*leng^2 -22*leng   4*leng^2];
  
% lumped mass matrix

 elseif ipt==2

    m=zeros(4,4);
    mass=rho*area*leng;
    m=diag([mass/2  0  mass/2  0]);

% diagonal mass matrix

 else

    m=zeros(4,4);
    mass=rho*area*leng;
    m=mass*diag([1/2  leng^2/78  1/2  leng^2/78]);

end
 

