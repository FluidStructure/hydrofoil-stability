%----------------------------------------------------------------------------
% EX5.11.2.m                                                              
% to solve the transient two-dimensional Laplace's equation             
%   u,t = u,xx + u,yy ,  0 < x < 5, 0 < y < 2                                                                        
% boundary conditions:
%   u(0,y,t) = 100, u(5,y,t) = 100, 
%   u,y(x,0,t) = 0, u,y(x,2,t) = 0
% initial condition:
%   u(x,y,0) = 0 over the domain
% using bilinear rectangular elements and forward difference method
%(see Fig. 5.11.1 for the finite element mesh except for
% 8 rectangular elements instead of 16 triangular elements)
%
% Variable descriptions                                                      
%   k = element matrix for time-independent term (u,xx + u,yy)                                           
%   m = element matrix for time-dependent term (u,t)
%   f = element vector
%   kk = system matrix of k                                            
%   mm = system matrix of m
%   ff = system vector                                                 
%   gcoord = coordinate values of each node
%   nodes = nodal connectivity of each element
%   index = a vector containing system dofs associated with each element     
%   bcdof = a vector containing dofs associated with boundary conditions     
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'                                              
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
ftime=10;                % termination time
ntime=fix((ftime-stime)/deltt); % number of time increment

%---------------------------------------------
%  input data for nodal coordinate values
%  gcoord(i,j) where i->node no. and j->x or y
%---------------------------------------------

gcoord(1,1)=0.0;   gcoord(1,2)=0.0;   gcoord(2,1)=1.25;   gcoord(2,2)=0.0;
gcoord(3,1)=2.5;   gcoord(3,2)=0.0;   gcoord(4,1)=3.75;   gcoord(4,2)=0.0;
gcoord(5,1)=5.0;   gcoord(5,2)=0.0;   gcoord(6,1)=0.0;    gcoord(6,2)=1.0;
gcoord(7,1)=1.25;  gcoord(7,2)=1.0;   gcoord(8,1)=2.5;    gcoord(8,2)=1.0;
gcoord(9,1)=3.75;  gcoord(9,2)=1.0;   gcoord(10,1)=5.0;   gcoord(10,2)=1.0;
gcoord(11,1)=0.0;  gcoord(11,2)=2.0;  gcoord(12,1)=1.25;  gcoord(12,2)=2.0;
gcoord(13,1)=2.5;  gcoord(13,2)=2.0;  gcoord(14,1)=3.75;  gcoord(14,2)=2.0;
gcoord(15,1)=5.0;  gcoord(15,2)=2.0;  

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
bcval(1)=100;             % whose described value is 0 
bcdof(2)=5;             % second node is constrained
bcval(2)=100;             % whose described value is 0
bcdof(3)=6;             % third node is constrained
bcval(3)=100;             % whose described value is 0 
bcdof(4)=10;             % 4th node is constrained
bcval(4)=100;             % whose described value is 0
bcdof(5)=11;             % 5th node is constrained
bcval(5)=100;             % whose described value is 0 
bcdof(6)=15;             % 6th node is constrained
bcval(6)=100;             % whose described value is 0

%-----------------------------------------
%  initialization of matrices and vectors
%-----------------------------------------

ff=zeros(sdof,1);       % initialization of system vector
fn=zeros(sdof,1);       % initialization of effective system vector
fsol=zeros(sdof,1);     % solution vector
sol=zeros(2,ntime+1);   % vector containing time history solution
kk=zeros(sdof,sdof);    % initialization of system matrix
mm=zeros(sdof,sdof);    % initialization of system matrix
index=zeros(nnel*ndof,1);  % initialization of index vector

%-----------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%-----------------------------------------------------------------

for iel=1:nel           % loop for the total number of elements

nd(1)=nodes(iel,1); % 1st connected node for (iel)-th element
nd(2)=nodes(iel,2); % 2nd connected node for (iel)-th element
nd(3)=nodes(iel,3); % 3rd connected node for (iel)-th element
nd(4)=nodes(iel,3); % 4th connected node for (iel)-th element
x1=gcoord(nd(1),1); y1=gcoord(nd(1),2);% coord values of 1st node
x2=gcoord(nd(2),1); y2=gcoord(nd(2),2);% coord values of 2nd node
x3=gcoord(nd(3),1); y3=gcoord(nd(3),2);% coord values of 3rd node
x4=gcoord(nd(4),1); y4=gcoord(nd(4),2);% coord values of 4th node
xleng=x2-x1;        % element size in x-axis
yleng=y4-y1;        % element size in y-axis

index=feeldof(nd,nnel,ndof);% extract system dofs associated with element

k=felp2dr4(xleng,yleng); % compute element matrix (time-independent term)
m=felpt2r4(xleng,yleng); % compute element matrix (time-dependent term)

kk=feasmbl1(kk,k,index);  % assemble element matrices 
mm=feasmbl1(mm,m,index);  % assemble element matrices 

end

%-----------------------------
%   loop for time integration 
%-----------------------------

for in=1:sdof
fsol(in)=0.0;    % initial condition  
end

sol(1,1)=fsol(8);  % store time history solution for node no. 8
sol(2,1)=fsol(9);  % store time history solution for node no. 9

for it=1:ntime     % start loop for time integration

fn=deltt*ff+(mm-deltt*kk)*fsol;  % compute effective column vector

[mm,fn]=feaplyc2(mm,fn,bcdof,bcval); % apply boundary condition

fsol=mm\fn;     % solve the matrix equation 

sol(1,it+1)=fsol(8); % store time history solution for node no. 8
sol(2,it+1)=fsol(9); % store time history solution for node no. 9

end


%-----------------------------------
% analytical solution at node 8
%-----------------------------------

pi=4*atan(1);
esol=zeros(1,ntime+1);
ii=0;
for ti=0:deltt:ntime*deltt;
ii=ii+1;
for i=1:2:100
esol(ii)=esol(ii)+(1/i)*exp(-i*i*pi*pi*ti/25)*sin(i*pi/2);
end
end
esol=100-(100*4/pi)*esol;

%------------------------------------
% plot the solution at nodes 8 and 9
%------------------------------------

time=0:deltt:ntime*deltt;
plot(time,sol(1,:),'*',time,esol,'-');
xlabel('Time')
ylabel('Solution at nodes')

%---------------------------------------------------------------

