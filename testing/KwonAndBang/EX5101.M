%----------------------------------------------------------------------------
% EX5.10.1.m                                                              
% to solve the axisymmetric Laplace equation given as            
%   u,rr + (u,r)/r + u,zz =0,  4 < r < 6, 0 < z < 1                                                                        
%   u(4,z) = 100, u,r(6,z) = 20 
%   u,z(r,0) = 0, u,z(r,1) = 0
% using linear triangular elements 
%(see Fig. 5.10.1 for the finite element mesh)
%
% Variable descriptions                                                      
%   k = element matrix                                             
%   f = element vector
%   kk = system matrix                                             
%   ff = system vector                                                 
%   gcoord = coordinate values of each node
%   nodes = nodal connectivity of each element
%   index = a vector containing system dofs associated with each element     
%   bcdof = a vector containing dofs associated with boundary conditions     
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'                                              
%----------------------------------------------------------------------------            

%------------------------------------
%  input data for control parameters
%------------------------------------

clear
nel=10;                  % number of elements
nnel=3;                  % number of nodes per element
ndof=1;                  % number of dofs per node
nnode=12;                % total number of nodes in system
sdof=nnode*ndof;         % total system dofs  

%---------------------------------------------
%  input data for nodal coordinate values
%  gcoord(i,j) where i->node no. and j->x or y
%---------------------------------------------

gcoord(1,1)=4.0;   gcoord(1,2)=0.0;   gcoord(2,1)=4.0;   gcoord(2,2)=1.0;
gcoord(3,1)=4.4;   gcoord(3,2)=0.0;   gcoord(4,1)=4.4;   gcoord(4,2)=1.0;
gcoord(5,1)=4.8;   gcoord(5,2)=0.0;   gcoord(6,1)=4.8;   gcoord(6,2)=1.0;
gcoord(7,1)=5.2;   gcoord(7,2)=0.0;   gcoord(8,1)=5.2;   gcoord(8,2)=1.0;
gcoord(9,1)=5.6;   gcoord(9,2)=0.0;   gcoord(10,1)=5.6;  gcoord(10,2)=1.0;
gcoord(11,1)=6.0;  gcoord(11,2)=0.0;  gcoord(12,1)=6.0;  gcoord(12,2)=1.0;

%---------------------------------------------------------
%  input data for nodal connectivity for each element
%  nodes(i,j) where i-> element no. and j-> connected nodes
%---------------------------------------------------------

nodes(1,1)=1;    nodes(1,2)=4;    nodes(1,3)=2;
nodes(2,1)=1;    nodes(2,2)=3;    nodes(2,3)=4; 
nodes(3,1)=3;    nodes(3,2)=6;    nodes(3,3)=4;
nodes(4,1)=3;    nodes(4,2)=5;    nodes(4,3)=6; 
nodes(5,1)=5;    nodes(5,2)=8;    nodes(5,3)=6;
nodes(6,1)=5;    nodes(6,2)=7;    nodes(6,3)=8; 
nodes(7,1)=7;    nodes(7,2)=10;   nodes(7,3)=8;
nodes(8,1)=7;    nodes(8,2)=9;    nodes(8,3)=10; 
nodes(9,1)=9;    nodes(9,2)=12;   nodes(9,3)=10;
nodes(10,1)=9;   nodes(10,2)=11;  nodes(10,3)=12; 

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------

bcdof(1)=1;             % first node is constrained
bcval(1)=100;           % whose described value is 100
bcdof(2)=2;             % second node is constrained
bcval(2)=100;           % whose described value is 100

%-----------------------------------------
%  initialization of matrices and vectors
%-----------------------------------------

ff=zeros(sdof,1);       % initialization of system force vector
kk=zeros(sdof,sdof);    % initialization of system matrix
index=zeros(nnel*ndof,1);  % initialization of index vector

pi=4*atan(1);              % define pi
ff(11)=120*pi; ff(12)=120*pi;  % flux boundary condition

%-----------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%-----------------------------------------------------------------

for iel=1:nel           % loop for the total number of elements

nd(1)=nodes(iel,1); % 1st connected node for (iel)-th element
nd(2)=nodes(iel,2); % 2nd connected node for (iel)-th element
nd(3)=nodes(iel,3); % 3rd connected node for (iel)-th element
r1=gcoord(nd(1),1); z1=gcoord(nd(1),2);% coord values of 1st node
r2=gcoord(nd(2),1); z2=gcoord(nd(2),2);% coord values of 2nd node
r3=gcoord(nd(3),1); z3=gcoord(nd(3),2);% coord values of 3rd node

index=feeldof(nd,nnel,ndof);% extract system dofs associated with element

k=felpaxt3(r1,z1,r2,z2,r3,z3); % compute element matrix

kk=feasmbl1(kk,k,index);  % assemble element matrices 

end

%-----------------------------
%   apply boundary conditions
%-----------------------------

[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);

%----------------------------
%  solve the matrix equation
%----------------------------

fsol=kk\ff;   

%---------------------
% analytical solution
%---------------------

for i=1:nnode
r=gcoord(i,1); z=gcoord(i,2);
esol(i)=100-6*20*log(4)+6*20*log(r);
end

%------------------------------------
% print both exact and fem solutions
%------------------------------------

num=1:1:sdof;
store=[num' fsol esol']


%---------------------------------------------------------------

