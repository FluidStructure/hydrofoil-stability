%----------------------------------------------------------------------------
% Example 9.9.2                                                              
%   axisymmetric analysis of a thick walled cylinder 
%   subjected to internal pressure using isoparametric 
%   four-node elements           
%   (see Fig. 9.9.2 for the finite element mesh)
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
%   point2 = matrix containing sampling points
%   weight2 = matrix containing weighting coefficients
%   bcdof = a vector containing dofs associated with boundary conditions     
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'                                              
%----------------------------------------------------------------------------            

%------------------------------------
%  input data for control parameters
%------------------------------------

clear
nel=5;                   % number of elements
nnel=4;                  % number of nodes per element
ndof=2;                  % number of dofs per node
nnode=12;                % total number of nodes in system
sdof=nnode*ndof;         % total system dofs  
edof=nnel*ndof;          % degrees of freedom per element
emodule=28.0e6;          % elastic modulus
poisson=0.25;            % Poisson's ratio
nglx=2; ngly=2;          % 2x2 Gauss-Legendre quadrature
nglxy=nglx*ngly;         % number of sampling points per element

%---------------------------------------------
%  input data for nodal coordinate values
%  gcoord(i,j) where i->node no. and j->x or y
%---------------------------------------------

gcoord=[10.  0.; 10.  1.; 11.  0.; 11.  1.; 12.  0.; 12.  1.; 
13.  0.; 13.  1.; 14.  0.; 14.  1.; 15.  0.; 15.  1.];

%---------------------------------------------------------
%  input data for nodal connectivity for each element
%  nodes(i,j) where i-> element no. and j-> connected nodes
%---------------------------------------------------------

nodes=[1 3 4 2; 3 5 6 4; 5 7 8 6; 7 9 10 8; 9 11 12 10];

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------

bcdof=[2 4 6 8 10 12 14 16 18 20 22 24];   % axial movement constrained
bcval=[0 0 0 0 0 0 0 0 0 0 0 0];        % whose described values are 0 

%-----------------------------------------
%  initialization of matrices and vectors
%-----------------------------------------

ff=zeros(sdof,1);       % system force vector
kk=zeros(sdof,sdof);    % system matrix
disp=zeros(sdof,1);     % system displacement vector
eldisp=zeros(edof,1);   % element displacement vector
stress=zeros(nglxy,4);  % matrix containing stress components
strain=zeros(nglxy,4);  % matrix containing strain components
index=zeros(edof,1);    % index vector
kinmtx=zeros(4,edof);   % kinematic matrix
matmtx=zeros(4,4);      % constitutive matrix

%----------------------------
%  force vector
%----------------------------

pi=4.0*atan(1.0);                          % pi=3.141592

ff(1)=2e3*2*pi*10;      % force applied at node 1 in x-axis
ff(3)=2e3*2*pi*10;      % force applied at node 2 in x-axis

%-----------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%-----------------------------------------------------------------

[point2,weight2]=feglqd2(nglx,ngly);       % sampling points & weights
matmtx=fematiso(3,emodule,poisson);        % compute constitutive matrix

for iel=1:nel           % loop for the total number of elements

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

[shape,dhdr,dhds]=feisoq4(x,y);    % compute shape functions and
				   % derivatives at sampling point

jacob2=fejacob2(nnel,dhdr,dhds,xcoord,ycoord);  % compute Jacobian

detjacob=det(jacob2);                 % determinant of Jacobian
invjacob=inv(jacob2);                 % inverse of Jacobian matrix

[dhdx,dhdy]=federiv2(nnel,dhdr,dhds,invjacob); % derivatives w.r.t.
					       % physical coordinate

xcenter=0;
for i=1:nnel                         % x-coordinate value
xcenter=xcenter+shape(i)*xcoord(i);  % of the integration point 
end

kinmtx=fekineax(nnel,dhdx,dhdy,shape,xcenter);   % kinematic matrix

%------------------------------
%  compute element matrix
%------------------------------

k=k+2*pi*xcenter*kinmtx'*matmtx*kinmtx*wtx*wty*detjacob; 
					% element matrix


end
end                                   % end of numerical integration loop

index=feeldof(nd,nnel,ndof);% extract system dofs for the element

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
end

%--------------------------------
%  numerical integration
%--------------------------------

intp=0;
for intx=1:nglx
x=point2(intx,1);                  % sampling point in x-axis
wtx=weight2(intx,1);               % weight in x-axis
for inty=1:ngly
y=point2(inty,2);                  % sampling point in y-axis
wty=weight2(inty,2) ;              % weight in y-axis
intp=intp+1;

[shape,dhdr,dhds]=feisoq4(x,y);    % compute shape functions and
				   % derivatives at sampling point

jacob2=fejacob2(nnel,dhdr,dhds,xcoord,ycoord);  % compute Jacobian

detjacob=det(jacob2);                 % determinant of Jacobian
invjacob=inv(jacob2);                 % inverse of Jacobian matrix

[dhdx,dhdy]=federiv2(nnel,dhdr,dhds,invjacob); % derivatives w.r.t.
					       % physical coordinate

xcenter=0;
for i=1:nnel                         % x-coordinate value
xcenter=xcenter+shape(i)*xcoord(i);  % of the integration point 
end

kinmtx=fekineax(nnel,dhdx,dhdy,shape,xcenter);   % kinematic matrix

index=feeldof(nd,nnel,ndof);% extract system dofs associated with element

%-------------------------------------------------------
%  extract element displacement vector
%-------------------------------------------------------

for i=1:edof
eldisp(i)=disp(index(i));
end

estrain=kinmtx*eldisp;             % compute strains
estress=matmtx*estrain;            % compute stresses

for i=1:4
strain(intp,i)=estrain(i);          % store for each element
stress(intp,i)=estress(i);          % store for each element          
end

end
end                                 % end of integration loop

for j=1:nglxy
stresses=[ielp stress(j,:)]                     % print stresses
end

end

%---------------------------------------------------------------

