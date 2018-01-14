%----------------------------------------------------------------------------
% EX5.12.1.m                                                              
% to solve the three-dimensional Laplace equation  
% for a pyramid shape of domain           
% using four-node tetrahedral elements.
% Bottom face has essential boundary condition and the 
% side faces are insulated. 
%(see Fig. 5.12.1 for the finite element mesh)
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
nel=4;                   % number of elements
nnel=4;                  % number of nodes per element
ndof=1;                  % number of dofs per node
nnode=6;                 % total number of nodes in system
sdof=nnode*ndof;         % total system dofs  

%---------------------------------------------
%  input data for nodal coordinate values
%  gcoord(i,j) where i->node no. and j->x or y
%---------------------------------------------

gcoord(1,1)=0.0;   gcoord(1,2)=0.0;   gcoord(1,3)=0.0;  
gcoord(2,1)=1.0;   gcoord(2,2)=0.0;   gcoord(2,3)=0.0; 
gcoord(3,1)=0.5;   gcoord(3,2)=0.5;   gcoord(3,3)=0.0;  
gcoord(4,1)=0.0;   gcoord(4,2)=1.0;   gcoord(4,3)=0.0;
gcoord(5,1)=1.0;   gcoord(5,2)=1.0;   gcoord(5,3)=0.0; 
gcoord(6,1)=0.5;   gcoord(6,2)=0.5;   gcoord(6,3)=1.0; 

%---------------------------------------------------------
%  input data for nodal connectivity for each element
%  nodes(i,j) where i-> element no. and j-> connected nodes
%---------------------------------------------------------

nodes(1,1)=4;  nodes(1,2)=1;  nodes(1,3)=3;  nodes(1,4)=6;
nodes(2,1)=1;  nodes(2,2)=2;  nodes(2,3)=3;  nodes(2,4)=6;
nodes(3,1)=2;  nodes(3,2)=5;  nodes(3,3)=3;  nodes(3,4)=6;
nodes(4,1)=5;  nodes(4,2)=4;  nodes(4,3)=3;  nodes(4,4)=6;

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------
 
bcdof(1)=1;             % 1st node is constrained
bcval(1)=0;             % whose described value is 0 
bcdof(2)=2;             % 2nd node is constrained
bcval(2)=20;            % whose described value is 20 
bcdof(3)=3;             % 3rd node is constrained
bcval(3)=150;           % whose described value is 150
bcdof(4)=4;             % 4th node is constrained
bcval(4)=100;           % whose described value is 100 
bcdof(5)=5;             % 5th node is constrained
bcval(5)=50;            % whose described value is 50 

%-----------------------------------------
%  initialization of matrices and vectors
%-----------------------------------------

ff=zeros(sdof,1);       % initialization of system force vector
kk=zeros(sdof,sdof);    % initialization of system matrix
index=zeros(nnel*ndof,1);  % initialization of index vector

%-----------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%-----------------------------------------------------------------

for iel=1:nel           % loop for the total number of elements

nd(1)=nodes(iel,1); % 1st connected node for (iel)-th element
nd(2)=nodes(iel,2); % 2nd connected node for (iel)-th element
nd(3)=nodes(iel,3); % 3rd connected node for (iel)-th element
nd(4)=nodes(iel,4); % 4th connected node for (iel)-th element
x(1)=gcoord(nd(1),1); y(1)=gcoord(nd(1),2);
z(1)=gcoord(nd(1),3);            % coord of 1st node
x(2)=gcoord(nd(2),1); y(2)=gcoord(nd(2),2);
z(2)=gcoord(nd(2),3);            % coord of 2nd node
x(3)=gcoord(nd(3),1); y(3)=gcoord(nd(3),2);
z(3)=gcoord(nd(3),3);            % coord of 3rd node
x(4)=gcoord(nd(4),1); y(4)=gcoord(nd(4),2);
z(4)=gcoord(nd(4),3);            % coord of 4th node

index=feeldof(nd,nnel,ndof);% extract system dofs associated with element

k=felp3dt4(x,y,z); % compute element matrix

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

%------------------------------------
% print both exact and fem solutions
%------------------------------------

num=1:1:sdof;
store=[num' fsol]


%---------------------------------------------------------------

