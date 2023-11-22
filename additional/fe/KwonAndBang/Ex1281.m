%----------------------------------------------------------------------------%
% Example 12.8.1                                                             %
% to find the critical buckling loads of of a simply supported               %
% beam using Hermitian beam elements                                         %           %
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
tleng=1;         % total length of the beam
leng=tleng/nel;  % uniform mesh (equal size of elements)

kk=zeros(sdof,sdof);  % initialization of system stiffness matrix
kkg=zeros(sdof,sdof);  % initialization of system geomtric matrix 
index=zeros(nel*ndof,1);  % initialization of index vector

bcdof(1)=1;		% deflection at node 1 is constrained
bcdof(2)=sdof-1;	% deflection at the last node is constrained

for iel=1:nel    % loop for the total number of elements

index=feeldof1(iel,nnel,ndof);  % extract system dofs associated with element

[k,m]=febeam1(el,xi,leng,0,0,1); % compute element stiffness & mass matrix

[kg,kef]=febeambk(leng,0); % compute geometric stiffness matrix

kk=feasmbl1(kk,k,index); % assemble element stiffness matrices into system matrix

kkg=feasmbl1(kkg,kg,index); % assemble geometric matrices into system matrix

end

[kk,kkg]=feaplycs(kk,kkg,bcdof); % apply constraints

fsol=eig(kk,kkg)   % solve the eigenvalue problem



