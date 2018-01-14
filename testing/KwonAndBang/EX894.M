%----------------------------------------------------------------------------%
% Example 8.9.4                                                              %
% to solve a static beam deflection problem using mixed beam elements        %
%                                                                            %
% Problem description                                                        %
%                                                                            %
%   Find the deflection of a simply supported beam whose length is           %
%   20 inches. The beam has also elastic modulus of 10x10e6 psi and          %
%   moment of inertia of cross-section 1/12 inch^4 with unit width.          %
%   It is subjected to a center load of 100 lb. Use 5 elements for           %
%   one half of the beam due to symmetry.                                    %
%   (see Fig. 8.9.1 for the element discretization)                          %
%                                                                            %
% Variable descriptions                                                      %
%   k = element stiffness matrix                                             %
%   kk = system stiffness matrix                                             %
%   ff = system force vector                                                 %
%   index = a vector containing system dofs associated with each element     %
%   bcdof = a vector containing dofs associated with boundary conditions     %
%   bcval = a vector containing boundary condition values associated with    %
%           the dofs in 'bcdof'                                              %
%----------------------------------------------------------------------------% 

clear
nel=5;           % number of elements
nnel=2;          % number of nodes per element
ndof=2;          % number of dofs per node
nnode=(nnel-1)*nel+1;   % total number of nodes in system
sdof=nnode*ndof; % total system dofs  

bcdof(1)=1;      % bending moment at node 1 is constrained
bcval(1)=0;      % whose described value is 0 
bcdof(2)=2;      % deflection at node 1 is constrained 
bcval(2)=0;      % whose described value is 0

ff=zeros(sdof,1);     % initialization of system force vector
kk=zeros(sdof,sdof);  % initialization of system matrix
index=zeros(nel*ndof,1);  % initialization of index vector
ff(12)=-50;       % because a half of the load is applied due to symmetry

for iel=1:nel    % loop for the total number of elements

index=feeldof1(iel,nnel,ndof);  % extract system dofs associated with element

k=febeam4(10^7,0.083333,2,0,1,1,1); % compute element stiffness matrix

kk=feasmbl1(kk,k,index); % assemble each element matrix into system matrix

end

[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);  % apply the boundary conditions

fsol=kk\ff;   % solve the matrix equation and print

% analytical solution

e=10^7;  l=20;  xi=1/12;  P=100;

for i = 1:nnode

x=(i-1)*2;
c=P/(48*e*xi);
k=(i-1)*ndof+1;
esol(k+1)=c*(3*l^2-4*x^2)*x;
esol(k)=-50*x;

end

% print both exact and fem solutions


num=1:1:sdof;
store=[num' fsol esol']


