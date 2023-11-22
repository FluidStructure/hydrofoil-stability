%-----------------------------------------------------------------------------
% Example 6.6.2                                                              
%   to solve the two-dimensional Laplace's equation given as            
%   u,xx + u,yy =0,  0 < x < 5, 0 < y < 10                                                                        
%   u(x,0) = 0, u(x,10) = 100sin(pi*x/10), 
%   u(0,y) = 0, u,x(5,y) = 0
%   using isoparametric four-node quadrilateral elements 
%   (see Fig. 5.9.2 for the finite element mesh)
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
%   point2 - integration (or sampling) points                                             
%   weight2 - weighting coefficients                                             
%   nglx - number of integration points along x-axis                                                
%   ngly - number of integration points along y-axis
%   xcoord - x coordinate values of nodes
%   ycoord - y coordinate values of nodes
%   jacob2 - jacobian matrix
%   shape - four-node quadrilateral shape functions
%   dhdr - derivatives of shape functions w.r.t. natural coord. r
%   dhds - derivatives of shape functions w.r.t. natural coord. s
%   dhdx - derivatives of shape functions w.r.t. physical coord. x
%   dhdy - derivatives of shape functions w.r.t. physical coord. y
%-----------------------------------------------------------------------------            
clear
%------------------------------------
%  input data for control parameters
%------------------------------------

nel=16;                  % number of elements
nnel=4;                  % number of nodes per element
ndof=1;                  % number of dofs per node
nnode=25;                % total number of nodes in system
nglx=2; ngly=2;          % use 2x2 integration rule
sdof=nnode*ndof;         % total system dofs  
edof=nnel*ndof;          % dofs per element

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

nodes(1,1)=1;    nodes(1,2)=2;    nodes(1,3)=7;    nodes(1,4)=6;
nodes(2,1)=2;    nodes(2,2)=3;    nodes(2,3)=8;    nodes(2,4)=7;
nodes(3,1)=3;    nodes(3,2)=4;    nodes(3,3)=9;    nodes(3,4)=8;
nodes(4,1)=4;    nodes(4,2)=5;    nodes(4,3)=10;   nodes(4,4)=9;
nodes(5,1)=6;    nodes(5,2)=7;    nodes(5,3)=12;   nodes(5,4)=11;
nodes(6,1)=7;    nodes(6,2)=8;    nodes(6,3)=13;   nodes(6,4)=12;
nodes(7,1)=8;    nodes(7,2)=9;    nodes(7,3)=14;   nodes(7,4)=13;
nodes(8,1)=9;    nodes(8,2)=10;   nodes(8,3)=15;   nodes(8,4)=14;
nodes(9,1)=11;   nodes(9,2)=12;   nodes(9,3)=17;   nodes(9,4)=16;
nodes(10,1)=12;  nodes(10,2)=13;  nodes(10,3)=18;  nodes(10,4)=17;
nodes(11,1)=13;  nodes(11,2)=14;  nodes(11,3)=19;  nodes(11,4)=18;
nodes(12,1)=14;  nodes(12,2)=15;  nodes(12,3)=20;  nodes(12,4)=19;
nodes(13,1)=16;  nodes(13,2)=17;  nodes(13,3)=22;  nodes(13,4)=21;
nodes(14,1)=17;  nodes(14,2)=18;  nodes(14,3)=23;  nodes(14,4)=22;
nodes(15,1)=18;  nodes(15,2)=19;  nodes(15,3)=24;  nodes(15,4)=23;
nodes(16,1)=19;  nodes(16,2)=20;  nodes(16,3)=25;  nodes(16,4)=24;

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
bcdof(10)=22;           % 22nd node is constrained
bcval(10)=38.2683;      % whose described value is 38.2683
bcdof(11)=23;           % 23rd node is constrained
bcval(11)=70.7107;      % whose described value is 70.7107
bcdof(12)=24;           % 24th node is constrained
bcval(12)=92.3880;      % whose described value is 92.3880
bcdof(13)=25;           % 25th node is constrained
bcval(13)=100;          % whose described value is 100

%-----------------------------------------
%  initialization of matrices and vectors
%-----------------------------------------

ff=zeros(sdof,1);       % initialization of system force vector
kk=zeros(sdof,sdof);    % initialization of system matrix
index=zeros(nnel*ndof,1);  % initialization of index vector

%-----------------------------------------------------------
%  loop for computation and assembly of element matrices
%-----------------------------------------------------------

[point2,weight2]=feglqd2(nglx,ngly);  % sampling points & weights

for iel=1:nel               % loop for the total number of elements

for i=1:nnel
nd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
xcoord(i)=gcoord(nd(i),1);  % extract x value of the node
ycoord(i)=gcoord(nd(i),2);  % extract y value of the node
end

k=zeros(edof,edof);         % initialization of element matrix to zero

%--------------------------------
%  numerical integration
%--------------------------------

for intx=1:nglx
x=point2(intx,1);                  % sampling point in x-axis
wtx=weight2(intx,1);               % weight in x-axis
for inty=1:ngly
y=point2(inty,2);                  % sampling point in y-axis
wty=weight2(inty,2) ;              % weight in y-axis

[shape,dhdr,dhds]=feisoq4(x,y); % compute shape functions and
                                      % derivatives at sampling point

jacob2=fejacob2(nnel,dhdr,dhds,xcoord,ycoord);  % compute Jacobian

detjacob=det(jacob2);                 % determinant of Jacobian
invjacob=inv(jacob2);                 % inverse of Jacobian matrix

[dhdx,dhdy]=federiv2(nnel,dhdr,dhds,invjacob); % derivatives w.r.t.
                                               % physical coordinate

%------------------------------
%  compute element matrix
%------------------------------

for i=1:edof
for j=1:edof
k(i,j)=k(i,j)+(dhdx(i)*dhdx(j)+dhdy(i)*dhdy(j))*wtx*wty*detjacob;
end
end

end
end                                   % end of numerical integration loop

index=feeldof(nd,nnel,ndof);% extract system dofs associated with element

%----------------------------------
% assemble element matrices 
%----------------------------------

kk=feasmbl1(kk,k,index);    

end                        % end of element loops

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

