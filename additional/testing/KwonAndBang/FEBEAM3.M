function [k,m]=febeam3(el,sh,leng,heig,width,rho)

%--------------------------------------------------------------
%  Purpose:
%     Stiffness and mass matrices for beam element with displacement
%     degrees of freedom only
%     nodal dof {u_1^b u_1^t v_1 u_2^b u_2^t v_2}
%
%  Synopsis:
%     [k,m]=febeam1(el,sh,leng,heig,rho,area,ipt) 
%
%  Variable Description:
%     k - element stiffness matrix (size of 6x6)    
%     m - element mass matrix (size of 6x6)
%     el - elastic modulus 
%     sh - shear modulus
%     leng - element length
%     heig - element thickness
%     width - width of the beam element
%     rho - mass density of the beam element (mass per unit volume)
%           lumped mass matrix only
%---------------------------------------------------------------

% stiffness matrix

 a1=(sh*leng*width)/(4*heig);
 a2=(sh*heig*width)/leng;
 a3=(el*heig*width)/(6*leng);
 a4=sh*width/2;
 k= [ a1+2*a3   -a1+a3     a4   a1-2*a3   -a1-a3    -a4;...
     -a1+a3      a1+2*a3  -a4  -a1-a3      a1-2*a3   a4;...
      a4        -a4        a2   a4        -a4       -a2;...
      a1-2*a3   -a1-a3     a4   a1+2*a3   -a1+a3    -a4;...
     -a1-a3      a1-2*a3  -a4  -a1+a3      a1+2*a3     a4;...
     -a4         a4       -a2  -a4         a4        a2];
   

% lumped mass matrix

    m=zeros(6,6);
    mass=rho*heig*width*leng/4;
    m=mass*diag([1  1   2   1   1   2]);


 

