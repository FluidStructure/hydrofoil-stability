%----------------------------------------------------------------------------%
% Example 8.10.1                                                             %
% to find the natural frequencies and mode shapes of a beam using            %
% Hermitian beam elements                                                    %
%                                                                            %
% Problem description                                                        %
%   Find the natural frequencies and mode shapes of a free beam of length    %
%   1. It has a cross-section 1 by 1 and it has also mass density of 1.      %
%   The elastic modulus of the beam is 12.                                   %
%   Use 4 elements to model the whole beam such that nonsymmetric            %
%   mode shapes can be included. use also consistent mass matrices.          %                                                                 
%                                                                            %
% Variable descriptions                                                      %
%   k = element stiffness matrix                                             %
%   m = element mass matrix                                                  %   
%   kk = system stiffness matrix                                             %
%   mm = system mass matrix                                                  %
%   index = a vector containing system dofs associated with each element     %
%   bcdof = a vector containing dofs associated with boundary conditions     %
%   bcval = a vector containing boundary condition values associated with    %
%           the dofs in 'bcdof'                                              %
%----------------------------------------------------------------------------%            

clear
nel=4;           % number of elements
nnel=2;          % number of nodes per element
ndof=2;          % number of dofs per node
nnode=(nnel-1)*nel+1;   % total number of nodes in system
sdof=nnode*ndof; % total system dofs  

el=12;           % elastic modulus
xi=1/12;         % moment of inertia of cross-section
rho=1;           % mass density
tleng=1;         % total length of the beam
leng=tleng/nel;  % uniform mesh (equal size of elements)
area=1;          % cross-sectional area

kk=zeros(sdof,sdof);  % initialization of system stiffness matrix
mm=zeros(sdof,sdof);  % initialization of system mass matrix 
index=zeros(nel*ndof,1);  % initialization of index vector


for iel=1:nel    % loop for the total number of elements

index=feeldof1(iel,nnel,ndof);  % extract system dofs associated with element

[k,m]=febeam1(el,xi,leng,area,rho,1); % compute element stiffness & mass matrix

kk=feasmbl1(kk,k,index); % assemble element stiffness matrices into system matrix

mm=feasmbl1(mm,m,index); % assemble element mass matrices into system matrix

end


fsol=eig(kk,mm);   % solve the eigenvalue problem
fsol=sqrt(fsol)

%_____________________________________________________________________
