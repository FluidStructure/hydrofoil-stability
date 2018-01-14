function [kg,kef]=febeambk(leng,elfd)

%-------------------------------------------------------------------
%  Purpose:
%     Geometric stiffness matrix and stiffness matrix for 
%     elastic foundation using Hermitian beam element
%     nodal dof {v_1 theta_1 v_2 theta_2}
%
%  Synopsis:
%     [kg,kef]=febeambk(el,xi,leng,area,rho,ipt) 
%
%  Variable Description:
%     kg - element geometric stiffness matrix (size of 4x4)   
%     kef - stiffness matrix for elastic foundation(size of 4x4)
%     leng - element length
%     elfd - spring constant of elastic foundation (force/length^2)
%-------------------------------------------------------------------

% geometric stiffness matrix

 c1=1/(30*leng);
 kg=c1*[36    3*leng   -36       3*leng;...
      3*leng  4*leng^2 -3*leng  -leng^2;...
     -36     -3*leng    36      -3*leng;...
      3*leng -leng^2   -3*leng   4*leng^2];
   
% stiffness matrix for elastic foundation

 d1=elfd*leng/420;
 kef=d1*[156      22*leng   54       -13*leng;...
          22*leng  4*leng^2  13*leng  -3*leng^2;...
          54       13*leng   156      -22*leng;...
         -13*leng -3*leng^2 -22*leng   4*leng^2];
  


