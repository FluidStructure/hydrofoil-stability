%----------------------------------------------------------------------------
% EX5.9.1.m                                                              
% to solve the two-dimensional Laplace's equation given as            
%   u,xx + u,yy =0,  0 < x < 5, 0 < y < 10                                                                        
%   u(x,0) = 0, u(x,10) = 100sin(pi*x/10), 
%   u(0,y) = 0, u,x(5,y) = 0
% using linear triangular elements 
%(see Fig. 5.9.1 for the finite element mesh)
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
nel=32;                  % number of elements
nnel=3;                  % number of nodes per element
ndof=1;                  % number of dofs per node
nnode=25;                % total number of nodes in system
sdof=nnode*ndof;         % total system dofs  

%---------------------------------------------
%  input data for nodal coordinate values
%  gcoord(i,j) where i->node no. and j->x or y
%---------------------------------------------

gcoord(1,1)=0.0;   gcoord(1,2)=0.0;   gcoord(2,1)=1.25;   gcoord(2,2)=0.0;
gcoord(3,1)=2.5;   gcoord(3,2)=0.0;   gcoord(4,1)=3.75;   gcoord(4,2)=0.0;
gcoord(5,1)=5.0;   gcoord(5,2)=0.0;   gcoord(6,1)=0.0;    gcoord(6,2)=2.5;
gcoord(7,1)=1.25;  gcoord(7,2)=2.5;   gcoord(8,1)=2.5;    gcoord(8,2)=2.5;
gcoord(9,1)=3.75;  gcoord(9,2)=2.5;   gcoord(10,1)=5.0;   gcoord(10,2)=2.5;
gcoord(11,1)=0.0;  gcoord(11,2)=5.0;  gcoord(12,1)=1.25;  gcoord(12,2)=5.0;
gcoord(13,1)=2.5;  gcoord(13,2)=5.0;  gcoord(14,1)=3.75;  gcoord(14,2)=5.0;
gcoord(15,1)=5.0;  gcoord(15,2)=5.0;  gcoord(16,1)=0.0;   gcoord(16,2)=7.5;
gcoord(17,1)=1.25; gcoord(17,2)=7.5;  gcoord(18,1)=2.5;   gcoord(18,2)=7.5;
gcoord(19,1)=3.75; gcoord(19,2)=7.5;  gcoord(20,1)=5.0;   gcoord(20,2)=7.5;
gcoord(21,1)=0.0;  gcoord(21,2)=10.;  gcoord(22,1)=1.25;  gcoord(22,2)=10.;
gcoord(23,1)=2.5;  gcoord(23,2)=10.;  gcoord(24,1)=3.75;  gcoord(24,2)=10.;
gcoord(25,1)=5.0;  gcoord(25,2)=10.;  

%---------------------------------------------------------
%  input data for nodal connectivity for each element
%  nodes(i,j) where i-> element no. and j-> connected nodes
%---------------------------------------------------------

nodes(1,1)=1;    nodes(1,2)=2;    nodes(1,3)=7;
nodes(2,1)=2;    nodes(2,2)=3;    nodes(2,3)=8; 
nodes(3,1)=3;    nodes(3,2)=4;    nodes(3,3)=9;
nodes(4,1)=4;    nodes(4,2)=5;    nodes(4,3)=10; 
nodes(5,1)=1;    nodes(5,2)=7;    nodes(5,3)=6;
nodes(6,1)=2;    nodes(6,2)=8;    nodes(6,3)=7; 
nodes(7,1)=3;    nodes(7,2)=9;    nodes(7,3)=8;
nodes(8,1)=4;    nodes(8,2)=10;   nodes(8,3)=9; 
nodes(9,1)=6;    nodes(9,2)=7;    nodes(9,3)=12;
nodes(10,1)=7;   nodes(10,2)=8;   nodes(10,3)=13; 
nodes(11,1)=8;   nodes(11,2)=9;   nodes(11,3)=14;
nodes(12,1)=9;   nodes(12,2)=10;  nodes(12,3)=15; 
nodes(13,1)=6;   nodes(13,2)=12;  nodes(13,3)=11;
nodes(14,1)=7;   nodes(14,2)=13;  nodes(14,3)=12; 
nodes(15,1)=8;   nodes(15,2)=14;  nodes(15,3)=13;
nodes(16,1)=9;   nodes(16,2)=15;  nodes(16,3)=14; 
nodes(17,1)=11;  nodes(17,2)=12;  nodes(17,3)=17;
nodes(18,1)=12;  nodes(18,2)=13;  nodes(18,3)=18; 
nodes(19,1)=13;  nodes(19,2)=14;  nodes(19,3)=19;
nodes(20,1)=14;  nodes(20,2)=15;  nodes(20,3)=20; 
nodes(21,1)=11;  nodes(21,2)=17;  nodes(21,3)=16;
nodes(22,1)=12;  nodes(22,2)=18;  nodes(22,3)=17; 
nodes(23,1)=13;  nodes(23,2)=19;  nodes(23,3)=18;
nodes(24,1)=14;  nodes(24,2)=20;  nodes(24,3)=19; 
nodes(25,1)=16;  nodes(25,2)=17;  nodes(25,3)=22;
nodes(26,1)=17;  nodes(26,2)=18;  nodes(26,3)=23; 
nodes(27,1)=18;  nodes(27,2)=19;  nodes(27,3)=24;
nodes(28,1)=19;  nodes(28,2)=20;  nodes(28,3)=25; 
nodes(29,1)=16;  nodes(29,2)=22;  nodes(29,3)=21;
nodes(30,1)=17;  nodes(30,2)=23;  nodes(30,3)=22; 
nodes(31,1)=18;  nodes(31,2)=24;  nodes(31,3)=23;
nodes(32,1)=19;  nodes(32,2)=25;  nodes(32,3)=24; 

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------

bcdof(1)=1;             % first node is constrained
bcval(1)=0;             % whose described value is 0 
bcdof(2)=2;             % second node is constrained
bcval(2)=0;             % whose described value is 0
bcdof(3)=3;             % third node is constrained
bcval(3)=0;             % whose described value is 0 
bcdof(4)=4;             % 4th node is constrained
bcval(4)=0;             % whose described value is 0
bcdof(5)=5;             % 5th node is constrained
bcval(5)=0;             % whose described value is 0 
bcdof(6)=6;             % 6th node is constrained
bcval(6)=0;             % whose described value is 0
bcdof(7)=11;            % 11th node is constrained
bcval(7)=0;             % whose described value is 0 
bcdof(8)=16;            % 16th node is constrained
bcval(8)=0;             % whose described value is 0
bcdof(9)=21;            % 21st node is constrained
bcval(9)=0;             % whose described value is 0 
bcdof(10)=22;           % second node is constrained
bcval(10)=38.2683;      % whose described value is 38.2683
bcdof(11)=23;           % third node is constrained
bcval(11)=70.7107;      % whose described value is 70.7107
bcdof(12)=24;           % 4th node is constrained
bcval(12)=92.3880;      % whose described value is 92.3880
bcdof(13)=25;           % 5th node is constrained
bcval(13)=100;          % whose described value is 100

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
x1=gcoord(nd(1),1); y1=gcoord(nd(1),2);% coord values of 1st node
x2=gcoord(nd(2),1); y2=gcoord(nd(2),2);% coord values of 2nd node
x3=gcoord(nd(3),1); y3=gcoord(nd(3),2);% coord values of 3rd node

index=feeldof(nd,nnel,ndof);% extract system dofs associated with element

k=felp2dt3(x1,y1,x2,y2,x3,y3); % compute element matrix

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
x=gcoord(i,1); y=gcoord(i,2);
esol(i)=100*sinh(0.31415927*y)*sin(0.31415927*x)/sinh(3.1415927); 
end

%------------------------------------
% print both exact and fem solutions
%------------------------------------

num=1:1:sdof;
store=[num' fsol esol']


%---------------------------------------------------------------

