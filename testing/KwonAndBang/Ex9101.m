%----------------------------------------------------------------------------
% Example 9.10.1                                                              
%   three-dimensional analysis of a cube using isoparametric 
%   eight-node elements           
%   (see Fig. 9.10.1 for the finite element mesh)
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
%   point3 = matrix containing sampling points
%   weight3 = matrix containing weighting coefficients
%   bcdof = a vector containing dofs associated with boundary conditions     
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'                                              
%----------------------------------------------------------------------------            

%------------------------------------
%  input data for control parameters
%------------------------------------

clear
nel=1;                   % number of elements
nnel=8;                  % number of nodes per element
ndof=3;                  % number of dofs per node
nnode=8;                 % total number of nodes in system
sdof=nnode*ndof;         % total system dofs  
edof=nnel*ndof;          % degrees of freedom per element
emodule=1e5;             % elastic modulus
poisson=0.3;             % Poisson's ratio
nglx=2; ngly=2; nglz=2;  % 2x2x2 Gauss-Legendre quadrature
nglxy=nglx*ngly*nglz;    % number of sampling points per element

%---------------------------------------------
%  input data for nodal coordinate values
%  gcoord(i,j) where i->node no. and j->x or y
%---------------------------------------------

gcoord=[0.0  0.0  0.0; 1.0  0.0  0.0; 1.0  1.0  0.0; 0.0  1.0  0.0;  
0.0  0.0  1.0; 1.0  0.0  1.0; 1.0  1.0  1.0; 0.0  1.0  1.0];

%---------------------------------------------------------
%  input data for nodal connectivity for each element
%  nodes(i,j) where i-> element no. and j-> connected nodes
%---------------------------------------------------------

nodes=[1 2 3 4 5 6 7 8];

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------

bcdof=[1 2 3 5 6 9 12];        % first four dofs are constrained
bcval=[0 0 0 0 0 0 0];        % whose described values are 0 

%-----------------------------------------
%  initialization of matrices and vectors
%-----------------------------------------

ff=zeros(sdof,1);       % system force vector
kk=zeros(sdof,sdof);    % system matrix
disp=zeros(sdof,1);     % system displacement vector
eldisp=zeros(edof,1);   % element displacement vector
stress=zeros(nglxy,6);  % matrix containing stress components
strain=zeros(nglxy,6);  % matrix containing strain components
index=zeros(edof,1);    % index vector
kinmtx=zeros(6,edof);   % kinematic matrix
matmtx=zeros(6,6);      % constitutive matrix

%----------------------------
%  force vector
%----------------------------

ff(15)=250;              % force applied at node 5 in z-axis
ff(18)=250;              % force applied at node 6 in z-axis
ff(21)=250;              % force applied at node 7 in z-axis
ff(24)=250;              % force applied at node 8 in z-axis

%-----------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%-----------------------------------------------------------------

[point3,weight3]=feglqd3(nglx,ngly,ngly);  % sampling points & weights
matmtx=fematiso(4,emodule,poisson);        % compute constitutive matrix

for iel=1:nel           % loop for the total number of elements

for i=1:nnel
nd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
xcoord(i)=gcoord(nd(i),1);  % extract x value of the node
ycoord(i)=gcoord(nd(i),2);  % extract y value of the node
zcoord(i)=gcoord(nd(i),3);  % extract z value of the node
end

k=zeros(edof,edof);         % initialization of element matrix to zero

%--------------------------------
%  numerical integration
%--------------------------------

for intx=1:nglx
x=point3(intx,1);                  % sampling point in x-axis
wtx=weight3(intx,1);               % weight in x-axis
for inty=1:ngly
y=point3(inty,2);                  % sampling point in y-axis
wty=weight3(inty,2) ;              % weight in y-axis
for intz=1:nglz
z=point3(intz,3);                  % sampling point in z-axis
wtz=weight3(intz,3) ;              % weight in z-axis

[shape,dhdr,dhds,dhdt]=feisos8(x,y,z);  % compute shape functions 
                              % and derivatives at sampling point

jacob3=fejacob3(nnel,dhdr,dhds,dhdt,xcoord,ycoord,zcoord);  
                                               % compute Jacobian

detjacob=det(jacob3);                 % determinant of Jacobian
invjacob=inv(jacob3);                 % inverse of Jacobian matrix

[dhdx,dhdy,dhdz]=federiv3(nnel,dhdr,dhds,dhdt,invjacob); 
                             % derivatives w.r.t. physical coordinate

kinmtx=fekine3d(nnel,dhdx,dhdy,dhdz);          % compute kinematic matrix

%------------------------------
%  compute element matrix
%------------------------------

k=k+kinmtx'*matmtx*kinmtx*wtx*wty*wtz*detjacob;    % element matrix

end
end
end                                   % end of numerical integration loop

index=feeldof(nd,nnel,ndof);% extract system dofs associated with element

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

num=1:1:sdof;
displace=[num' disp]                  % print nodal displacements

%---------------------------------------
%  element stress computation
%---------------------------------------

for ielp=1:nel           % loop for the total number of elements

for i=1:nnel
nd(i)=nodes(ielp,i);        % extract connected node for (iel)-th element
xcoord(i)=gcoord(nd(i),1);  % extract x value of the node
ycoord(i)=gcoord(nd(i),2);  % extract y value of the node
zcoord(i)=gcoord(nd(i),3);  % extract z value of the node
end

%--------------------------------
%  numerical integration
%--------------------------------
intp=0;
for intx=1:nglx
x=point3(intx,1);                  % sampling point in x-axis
wtx=weight3(intx,1);               % weight in x-axis
for inty=1:ngly
y=point3(inty,2);                  % sampling point in y-axis
wty=weight3(inty,2) ;              % weight in y-axis
for intz=1:nglz
z=point3(intz,3);                  % sampling point in z-axis
wtz=weight3(intz,3) ;              % weight in z-axis
intp=intp+1;

[shape,dhdr,dhds,dhdt]=feisos8(x,y,z);  % compute shape functions 
                              % and derivatives at sampling point

jacob3=fejacob3(nnel,dhdr,dhds,dhdt,xcoord,ycoord,zcoord);  
                                               % compute Jacobian

detjacob=det(jacob3);                 % determinant of Jacobian
invjacob=inv(jacob3);                 % inverse of Jacobian matrix

[dhdx,dhdy,dhdz]=federiv3(nnel,dhdr,dhds,dhdt,invjacob); 
                             % derivatives w.r.t. physical coordinate

kinmtx=fekine3d(nnel,dhdx,dhdy,dhdz);          % compute kinematic matrix

index=feeldof(nd,nnel,ndof);% extract system dofs for the element

%-------------------------------------------------------
%  extract element displacement vector
%-------------------------------------------------------

for i=1:edof
eldisp(i)=disp(index(i));
end

estrain=kinmtx*eldisp;              % compute strains
estress=matmtx*estrain;             % compute stresses

for i=1:6
strain(intp,i)=estrain(i);          % store for each element
stress(intp,i)=estress(i);          % store for each element          
end


location=[ielp,intx,inty,intz]         % print location for stress
stress(intp,:)                    % print stress values

end
end
end                                 % end of integration loop


end

%---------------------------------------------------------------

