%----------------------------------------------------------------------------
% Example 9.8.1                                                              
%   plane stress analysis of a solid using linear triangular elements           
%   (see Fig. 9.8.1 for the finite element mesh)
%
% Variable descriptions                                                      
%   k = element matrix                                             
%   f = element vector
%   kk = system matrix                                             
%   ff = system vector                                                 
%   disp = system nodal displacement vector
%   eldisp = element nodal displacement vector
%   stress = matrix containing stresses
%   strain = matrix containing strains
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
nel=8;                   % number of elements
nnel=3;                  % number of nodes per element
ndof=2;                  % number of dofs per node
nnode=10;                % total number of nodes in system
sdof=nnode*ndof;         % total system dofs  
edof=nnel*ndof;          % degrees of freedom per element
emodule=100000.0;        % elastic modulus
poisson=0.3;             % Poisson's ratio

%---------------------------------------------
%  input data for nodal coordinate values
%  gcoord(i,j) where i->node no. and j->x or y
%---------------------------------------------

gcoord=[0.0  0.0; 0.0  1.0; 1.0  0.0; 1.0  1.0; 2.0  0.0; 
2.0  1.0; 3.0  0.0; 3.0  1.0; 4.0  0.0; 4.0  1.0];

%---------------------------------------------------------
%  input data for nodal connectivity for each element
%  nodes(i,j) where i-> element no. and j-> connected nodes
%---------------------------------------------------------

nodes=[1 3 4; 1 4 2; 3 5 6; 3 6 4; 
5 7 8; 5 8 6;7 9 10; 7 10 8];

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------

bcdof=[1 2 3];        % first three dofs are constrained
bcval=[0 0 0];        % whose described values are 0 

%-----------------------------------------
%  initialization of matrices and vectors
%-----------------------------------------

ff=zeros(sdof,1);       % system force vector
kk=zeros(sdof,sdof);    % system matrix
disp=zeros(sdof,1);     % system displacement vector
eldisp=zeros(edof,1);   % element displacement vector
stress=zeros(nel,3);    % matrix containing stress components
strain=zeros(nel,3);    % matrix containing strain components
index=zeros(edof,1);    % index vector
kinmtx2=zeros(3,edof);   % kinematic matrix
matmtx=zeros(3,3);      % constitutive matrix

%----------------------------
%  force vector
%----------------------------

ff(17)=500;              % force applied at node 9 in x-axis
ff(19)=500;              % force applied at node 10 in x-axis

%-----------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%-----------------------------------------------------------------

matmtx=fematiso(1,emodule,poisson);        % compute constitutive matrix

for iel=1:nel           % loop for the total number of elements

nd(1)=nodes(iel,1); % 1st connected node for (iel)-th element
nd(2)=nodes(iel,2); % 2nd connected node for (iel)-th element
nd(3)=nodes(iel,3); % 3rd connected node for (iel)-th element

x1=gcoord(nd(1),1); y1=gcoord(nd(1),2);% coord values of 1st node
x2=gcoord(nd(2),1); y2=gcoord(nd(2),2);% coord values of 2nd node
x3=gcoord(nd(3),1); y3=gcoord(nd(3),2);% coord values of 3rd node

index=feeldof(nd,nnel,ndof);% extract system dofs associated with element

%-------------------------------------------------------
%  find the derivatives of shape functions
%-------------------------------------------------------

area=0.5*(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2);  % area of triangule
area2=area*2;
dhdx=(1/area2)*[(y2-y3) (y3-y1) (y1-y2)];  % derivatives w.r.t. x-axis
dhdy=(1/area2)*[(x3-x2) (x1-x3) (x2-x1)];  % derivatives w.r.t. y-axis

kinmtx2=fekine2d(nnel,dhdx,dhdy);          % compute kinematic matrix

k=kinmtx2'*matmtx*kinmtx2*area;      % element stiffnes matrix

kk=feasmbl1(kk,k,index);  % assemble element matrices 

end

%-----------------------------
%   apply boundary conditions
%-----------------------------

[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);

%----------------------------
%  solve the matrix equation
%----------------------------

disp=kk\ff;   

%---------------------------------------
%  element stress computation
%---------------------------------------

for ielp=1:nel           % loop for the total number of elements

nd(1)=nodes(ielp,1); % 1st connected node for (iel)-th element
nd(2)=nodes(ielp,2); % 2nd connected node for (iel)-th element
nd(3)=nodes(ielp,3); % 3rd connected node for (iel)-th element

x1=gcoord(nd(1),1); y1=gcoord(nd(1),2);% coord values of 1st node
x2=gcoord(nd(2),1); y2=gcoord(nd(2),2);% coord values of 2nd node
x3=gcoord(nd(3),1); y3=gcoord(nd(3),2);% coord values of 3rd node

index=feeldof(nd,nnel,ndof);% extract system dofs associated with element

%-------------------------------------------------------
%  extract element displacement vector
%-------------------------------------------------------

for i=1:edof
eldisp(i)=disp(index(i));
end

area=0.5*(x1*y2+x2*y3+x3*y1-x1*y3-x2*y1-x3*y2);  % area of triangule
area2=area*2;
dhdx=(1/area2)*[(y2-y3) (y3-y1) (y1-y2)];  % derivatives w.r.t. x-axis
dhdy=(1/area2)*[(x3-x2) (x1-x3) (x2-x1)];  % derivatives w.r.t. y-axis

kinmtx2=fekine2d(nnel,dhdx,dhdy);          % compute kinematic matrix

estrain=kinmtx2*eldisp;             % compute strains
estress=matmtx*estrain;             % compute stresses

for i=1:3
strain(ielp,i)=estrain(i);          % store for each element
stress(ielp,i)=estress(i);          % store for each element          
end

end

%------------------------------------
% print fem solutions
%------------------------------------

num=1:1:sdof;
displace=[num' disp]                         % print nodal displacements

for i=1:nel
stresses=[i stress(i,:)]                     % print stresses
end

%---------------------------------------------------------------

