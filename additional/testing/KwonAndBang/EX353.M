%----------------------------------------------------------------------------
% EX3.5.3                                                              
% to solve the ordinary differential equation given as            
%   x^2 u'' - 2x u' - 4u = x^2,  10 < x < 20                                                                        
%   u(10) = 0  and  u(20) = 100
% using 10 linear elements
%
% Variable descriptions                                                      
%   k = element matrix                                             
%   f = element vector
%   kk = system matrix                                             
%   ff = system vector                                                 
%   index = a vector containing system dofs associated with each element     
%   bcdof = a vector containing dofs associated with boundary conditions     
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'                                              
%----------------------------------------------------------------------------            

%------------------------------------
%  input data for control parameters
%------------------------------------

clear
nel=10;                 % number of elements
nnel=2;                 % number of nodes per element
ndof=1;                 % number of dofs per node
nnode=11;               % total number of nodes in system
sdof=nnode*ndof;        % total system dofs  

%-----------------------------------------
%  input data for nodal coordinate values
%-----------------------------------------

gcoord(1)=10;  gcoord(2)=11;  gcoord(3)=12;  gcoord(4)=13;
gcoord(5)=14;  gcoord(6)=15;  gcoord(7)=16;  gcoord(8)=17;
gcoord(9)=18;  gcoord(10)=19;  gcoord(11)=20; 

%-----------------------------------------------------
%  input data for nodal connectivity for each element
%-----------------------------------------------------

nodes(1,1)=1;  nodes(1,2)=2;  nodes(2,1)=2;  nodes(2,2)=3;
nodes(3,1)=3;  nodes(3,2)=4;  nodes(4,1)=4;  nodes(4,2)=5;
nodes(5,1)=5;  nodes(5,2)=6;  nodes(6,1)=6;  nodes(6,2)=7;  
nodes(7,1)=7;  nodes(7,2)=8;  nodes(8,1)=8;  nodes(8,2)=9;
nodes(9,1)=9;  nodes(9,2)=10;  nodes(10,1)=10;  nodes(10,2)=11;

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------

bcdof(1)=1;             % first node is constrained
bcval(1)=0;             % whose described value is 0 
bcdof(2)=11;            % 11th node is constrained
bcval(2)=100;           % whose described value is 100 

%-----------------------------------------
%  initialization of matrices and vectors
%-----------------------------------------

ff=zeros(sdof,1);       % initialization of system force vector
kk=zeros(sdof,sdof);    % initialization of system matrix
index=zeros(nnel*ndof,1);  % initialization of index vector

%-----------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%-----------------------------------------------------------------

for iel=1:nel                % loop for the total number of elements

nl=nodes(iel,1); nr=nodes(iel,2); % extract nodes for (iel)-th element
xl=gcoord(nl); xr=gcoord(nr);% extract nodal coord values for the element
eleng=xr-xl;                 % element length
index=feeldof1(iel,nnel,ndof);% extract system dofs associated with element

k=feodex2l(xl,xr);                  % compute element matrix
f=fefx2l(xl,xr);                    % compute element vector
[kk,ff]=feasmbl2(kk,ff,k,f,index);  % assemble element matrices and vectors

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

esol(1)=0.0;
for i=2:nnode
x=gcoord(i);
esol(i)=0.00102*x^4-0.16667*x^2+64.5187/x; 
end

%------------------------------------
% print both exact and fem solutions
%------------------------------------

num=1:1:sdof;
store=[num' fsol esol']

%---------------------------------------------------------------------

