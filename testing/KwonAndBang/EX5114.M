%----------------------------------------------------------------------------
% EX5.11.4.m                                                              
% to solve the transient two-dimensional Laplace's equation             
%   a u,t = u,xx + u,yy ,  0 < x < 0.02, 0 < y < 0.01                                                                        
% boundary conditions:
%   u(0,y,t) = 300, u(0.02,y,t) = 300, 
%   u,y(x,0,t) = 0, u,y(x,0.01,t) = -(100/0.3)*(u-50)
% initial condition:
%   u(x,y,0) = 0 over the domain
% using bilinear rectangular elements and forward difference method
%(see Fig. 5.11.1 for the finite element mesh)
%
% Variable descriptions                                                      
%   k = element matrix for time-independent term (u,xx + u,yy)                                           
%   m = element matrix for time-dependent term (u,t)
%   f = element vector
%   kk = system matrix of k                                            
%   mm = system matrix of m
%   ff = system vector  
%   fn = effective system vector
%   fsol = solution vector
%   sol = time-history solution vector of selected nodes                                               
%   gcoord = coordinate values of each node
%   nodes = nodal connectivity of each element
%   index = a vector containing system dofs associated with each element     
%   bcdof = a vector containing dofs associated with boundary conditions     
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof' 
%   k1 = element matrix due to Cauchy-type flux
%   f1 = element vector due to flux boundary condition
%   index1 = index for nodal dofs with flux                                           
%----------------------------------------------------------------------------            
clear
%------------------------------------
%  input data for control parameters
%------------------------------------

nel=8;                   % number of elements
nnel=4;                  % number of nodes per element
ndof=1;                  % number of dofs per node
nnode=15;                % total number of nodes in system
sdof=nnode*ndof;         % total system dofs  
deltt=0.1;               % time step size for transient analysis
stime=0.0;               % initial time
ftime=1.0;               % termination time
ntime=fix((ftime-stime)/deltt); % number of time increment
a=4266.7;                % coefficient for the transient term
nf=4;                    % number of element boundaries with flux
nnels=2;                 % number of nodes per side of each element
%---------------------------------------------
%  input data for nodal coordinate values
%  gcoord(i,j) where i->node no. and j->x or y
%---------------------------------------------

gcoord(1,1)=0.0;     gcoord(1,2)=0.0;   
gcoord(2,1)=0.005;   gcoord(2,2)=0.0;
gcoord(3,1)=0.010;   gcoord(3,2)=0.0;   
gcoord(4,1)=0.015;   gcoord(4,2)=0.0;
gcoord(5,1)=0.020;   gcoord(5,2)=0.0;   
gcoord(6,1)=0.0;     gcoord(6,2)=0.005;
gcoord(7,1)=0.005;   gcoord(7,2)=0.005;   
gcoord(8,1)=0.010;   gcoord(8,2)=0.005;
gcoord(9,1)=0.015;   gcoord(9,2)=0.005;   
gcoord(10,1)=0.020;  gcoord(10,2)=0.005;
gcoord(11,1)=0.0;    gcoord(11,2)=0.01;  
gcoord(12,1)=0.005;  gcoord(12,2)=0.01;
gcoord(13,1)=0.010;  gcoord(13,2)=0.01;  
gcoord(14,1)=0.015;  gcoord(14,2)=0.01;
gcoord(15,1)=0.020;  gcoord(15,2)=0.01;  

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

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------

bcdof(1)=1;             % first node is constrained
bcval(1)=300;           % whose described value is 300 
bcdof(2)=5;             % second node is constrained
bcval(2)=300;           % whose described value is 300
bcdof(3)=6;             % third node is constrained
bcval(3)=300;           % whose described value is 300 
bcdof(4)=10;            % 4th node is constrained
bcval(4)=300;           % whose described value is 300
bcdof(5)=11;            % 5th node is constrained
bcval(5)=300;           % whose described value is 300 
bcdof(6)=15;            % 6th node is constrained
bcval(6)=300;           % whose described value is 300

%---------------------------------------------------------
%  input for flux boundary conditions
%  nflx(i,j) where i-> element no. and j-> two side nodes
%---------------------------------------------------------

nflx(1,1)=11;  nflx(1,2)=12;  % nodes on 1st element side with flux 
nflx(2,1)=12;  nflx(2,2)=13;  % nodes on 2nd element side with flux 
nflx(3,1)=13;  nflx(3,2)=14;  % nodes on 3rd element side with flux 
nflx(4,1)=14;  nflx(4,2)=15;  % nodes on 4th element side with flux 

b=333.3; c=50;  % Constants for Cauchy-type BC (du/dn=b(u-c))

%-----------------------------------------
%  initialization of matrices and vectors
%-----------------------------------------

ff=zeros(sdof,1);                % system vector
fn=zeros(sdof,1);                % effective system vector
fsol=zeros(sdof,1);              % solution vector
sol=zeros(1,ntime+1);            % time-history of a selected node
kk=zeros(sdof,sdof);             % of system matrix
mm=zeros(sdof,sdof);             % system matrix
index=zeros(nnel*ndof,1);        % index vector
f1=zeros(nnels*ndof,1);          % element flux vector
k1=zeros(nnels*ndof,nnels*ndof); % flux matrix
index1=zeros(nnels*ndof,1);      % flux index vector

%-----------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%-----------------------------------------------------------------

for iel=1:nel           % loop for the total number of elements

nd(1)=nodes(iel,1); % 1st connected node for (iel)-th element
nd(2)=nodes(iel,2); % 2nd connected node for (iel)-th element
nd(3)=nodes(iel,3); % 3rd connected node for (iel)-th element
nd(4)=nodes(iel,4); % 4th connected node for (iel)-th element
x1=gcoord(nd(1),1); y1=gcoord(nd(1),2);% coord values of 1st node
x2=gcoord(nd(2),1); y2=gcoord(nd(2),2);% coord values of 2nd node
x3=gcoord(nd(3),1); y3=gcoord(nd(3),2);% coord values of 3rd node
x4=gcoord(nd(4),1); y4=gcoord(nd(4),2);% coord values of 4th node
xleng=x2-x1;        % element size in x-axis
yleng=y4-y1;        % element size in y-axis

index=feeldof(nd,nnel,ndof);% extract system dofs associated with element

k=felp2dr4(xleng,yleng); % compute element matrix
m=a*felpt2r4(xleng,yleng); % compute element matrix

kk=feasmbl1(kk,k,index);  % assemble element matrices 
mm=feasmbl1(mm,m,index);  % assemble element matrices 

end

%--------------------------------------------------------
%  additional computation due to flux boundary condition
%--------------------------------------------------------

for ifx=1:nf

nds(1)=nflx(ifx,1);    % extract node with flux BC for (ifx)-th element
nds(2)=nflx(ifx,2);    % extarct node with flux BC for (ifx)-th element
x1=gcoord(nds(1),1); y1=gcoord(nds(1),2); % nodal coordinate
x2=gcoord(nds(2),1); y2=gcoord(nds(2),2); % nodal coordinate
eleng=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)); % element side length

index1=feeldof(nds,nnels,ndof); % find related system dofs

k1=b*feflxl2(eleng);  % compute element matrix due to flux
f1=b*c*fef1l(0,eleng); % compute element vector due to flux

[kk,ff]=feasmbl2(kk,ff,k1,f1,index1);  % assembly

end

%-----------------------------
%   loop for time integration 
%-----------------------------

for in=1:sdof
fsol(in)=300.0;    % initial condition
end

sol(1)=fsol(8);    % sol contains time-history solution at node 8

kk=mm+deltt*kk;

for it=1:ntime

fn=deltt*ff+mm*fsol;       % compute effective column vector

[kk,fn]=feaplyc2(kk,fn,bcdof,bcval); % apply boundary condition

fsol=kk\fn;     % solve the matrix equation 

sol(it+1)=fsol(8);  % sol contains time-history solution at node 8

end

%------------------------------------
% plot the solution at node 8
%------------------------------------

time=0:deltt:ntime*deltt;
plot(time,sol);
xlabel('Time')
ylabel('Solution at the center')

%---------------------------------------------------------------

