function [k,m]=fetruss1(el,leng,area,rho,ipt)

%--------------------------------------------------------------
%  Purpose:
%     Stiffness and mass matrices for the 1-d truss element
%     nodal dof {u_1 u_2}
%
%  Synopsis:
%     [k,m]=fetruss1(el,leng,area,rho,ipt) 
%
%  Variable Description:
%     k - element stiffness matrix (size of 4x4)   
%     m - element mass matrix (size of 4x4)
%     el - elastic modulus 
%     leng - element length
%     area - area of truss cross-section
%     rho - mass density (mass per unit volume)
%     ipt = 1 - consistent mass matrix
%         = 2 - lumped mass matrix
%--------------------------------------------------------------------------

% stiffness matrix 

 k= (area*el/leng)*[ 1  -1;...
                    -1   1];
   

% consistent mass matrix

 if ipt==1

    m=(rho*area*leng/6)*[ 2  1;...
                          1  2];
  
% lumped mass matrix

 else

    m=(rho*area*leng/2)*[ 1   0;...
                          0   1];

 end
 
