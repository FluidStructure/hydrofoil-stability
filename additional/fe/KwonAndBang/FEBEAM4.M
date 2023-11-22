function [k,m]=febeam4(el,xi,leng,sh,heig,rho,ipt)

%--------------------------------------------------------------
%  Purpose:
%     Stiffness and mass matrices for mixed beam element
%     bending moment and deflection as nodal degrees of freedom
%     nodal dof {M_1 v_1 M_2 v_2}
%
%  Synopsis:
%     [k,m]=febeam4(el,xi,leng,sh,heig,rho,ipt) 
%
%  Variable Description:
%     k - element stiffness matrix (size of 4x4)     
%     m - element mass matrix (size of 4x4)
%     el - elastic modulus 
%     xi - second moment of inertia of cross-section
%     leng - length of the beam element
%     sh - shear modulus              
%     heig - beam thickness
%     rho - mass density of the beam element (mass per unit volume)
%     ipt = 1           - lumped mass matrix
%         = otherwise   - diagonalized mass matrix
%---------------------------------------------------------------

% stiffness matrix

 if sh == 0

% thin beam (no shear deformation)

 k= [ leng/(3*el*xi)     1/leng   leng/(6*el*xi)    -1/leng;...
      1/leng             0          -1/leng          0;...       
      leng/(6*el*xi)    -1/leng   leng/(3*el*xi)     1/leng;...
     -1/leng             0           1/leng          0 ];

 else  

% thick beam (includes shear deformation)

 a=6/(5*sh*leng*heig);
 k= [ 1/(3*el*xi)+a   1/leng      1/(6*el*xi)-a  -1/leng;...
      1/leng          0          -1/leng          0;...       
      1/(6*el*xi)-a  -1/leng      1/(3*el*xi)+a   1/leng;...
     -1/leng          0           1/leng          0 ];

 end

% lumped mass matrix

if ipt==1

    m=zeros(4,4);
    mass=rho*heig*leng/2;
    m=diag([0  1  0  1]);

% diagonal mass matrix

else
        
    m=zeros(4,4);
    mass=rho*heig*leng/2;
    m=mass*diag([1  1  1  1]);

end
 

