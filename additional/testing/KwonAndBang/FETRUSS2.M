function [k,m]=fetruss2(el,leng,area,rho,beta,ipt)

%--------------------------------------------------------------
%  Purpose:
%     Stiffness and mass matrices for the 2-d truss element
%     nodal dof {u_1 v_1 u_2 v_2}
%
%  Synopsis:
%     [k,m]=fetruss2(el,leng,area,rho,beta,ipt) 
%
%  Variable Description:
%     k - element stiffness matrix (size of 4x4)   
%     m - element mass matrix (size of 4x4)
%     el - elastic modulus 
%     leng - element length
%     area - area of truss cross-section
%     rho - mass density (mass per unit volume)
%     beta - angle between the local and global axes                                 ipt = 1: consistent mass matrix
%            positive if the local axis is in the ccw direction from
%            the global axis
%     ipt = 1 - consistent mass matrix
%         = 2 - lumped mass matrix
%--------------------------------------------------------------------------

% stiffness matrix 

 c=cos(beta); s=sin(beta); 
 k= (area*el/leng)*[ c*c   c*s  -c*c  -c*s;...
                     c*s   s*s  -c*s  -s*s;...
                    -c*c  -c*s   c*c   c*s;...
                    -c*s  -s*s   c*s   s*s];
   

% consistent mass matrix

 if ipt==1

    m=(rho*area*leng/6)*[ 2*c*c+2*s*s  0  c*c+s*s  0;...
                          0  2*c*c+2*s*s  0  c*c+s*s;...
                          c*c+s*s  0  2*c*c+2*s*s  0;...
                          0  c*c+s*s  0  2*c*c+2*s*s];
  
% lumped mass matrix

 else

    m=(rho*area*leng/2)*[ c*c+s*s  0  0  0;...
                          0  c*c+s*s  0  0;...
                          0  0  c*c+s*s  0;...
                          0  0  0  c*c+s*s];

 end
 
