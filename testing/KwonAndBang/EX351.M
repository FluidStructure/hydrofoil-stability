%----------------------------------------------------------------------------
% EX3.5.1                                                              
% to solve the ordinary differential equation given as            
%   a u'' + b u' + c u = 1,  0 < x < 1                                                                        
%   u(0) = 0  and  u(1) = 0
% using 5 linear elements
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
nel=5;                  % number of elements
nnel=2;                 % number of nodes per element
ndof=1;                 % number of dofs per node
nnode=6;                % total number of nodes in system
sdof=nnode*ndof;        % total system dofs  

%-----------------------------------------
%  input data for nodal coordinate values
%-----------------------------------------

gcoord(1)=0.0;  gcoord(2)=0.2;  gcoord(3)=0.4;  gcoord(4)=0.6;
gcoord(5)=0.8;  gcoord(6)=1.0;   

%-----------------------------------------------------
%  input data for nodal connectivity for each element
%-----------------------------------------------------

nodes(1,1)=1;  nodes(1,2)=2;  nodes(2,1)=2;  nodes(2,2)=3;
nodes(3,1)=3;  nodes(3,2)=4;  nodes(4,1)=4;  nodes(4,2)=5;
nodes(5,1)=5;  nodes(5,2)=6;  

%-----------------------------------------
%  input data for coefficients of the ODE
%-----------------------------------------

acoef=1;                % coefficient 'a' of the diff eqn
bcoef=-3;               % coefficient 'b' of the diff eqn
ccoef=2;                % coefficient 'c' of the diff eqn    

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------

bcdof(1)=1;             % first node is constrained
bcval(1)=0;             % whose described value is 0 
bcdof(2)=6;             % 6th node is constrained
bcval(2)=0;             % whose described value is 0

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

nl=nodes(iel,1); nr=nodes(iel,2); % extract nodes for (iel)-th element
xl=gcoord(nl); xr=gcoord(nr);% extract nodal coord values for the element
eleng=xr-xl;            % element length
index=feeldof1(iel,nnel,ndof);% extract system dofs associated with element

k=feode2l(acoef,bcoef,ccoef,eleng); % compute element matrix
f=fef1l(xl,xr);                     % compute element vector
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

c1=0.5/exp(1);
c2=-0.5*(1+1/exp(1));
for i=1:nnode
x=gcoord(i);
esol(i)=c1*exp(2*x)+c2*exp(x)+1/2; 
end

%------------------------------------
% print both exact and fem solutions
%------------------------------------

num=1:1:sdof;
store=[num' fsol esol']


%---------------------------------------------------------------

