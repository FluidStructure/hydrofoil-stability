%----------------------------------------------------------------------------
% Example 10.9.1                                                              
%   A simply supported square plate is subjected to a concentrated load 
%   at the center. Find the deflection of the plate using 4 four-node
%   isoparametric elements of the shear deformable displacement
%   formulation. The size of the plate is 10 in. by 10 in. and its          
%   thickness is 0.1 in. It is made of steel and the applied 
%   force is 40 lb.
%   (see Fig. 10.9.1 for the finite element mesh)
%
% Variable descriptions                                                      
%   k = element matrix                                             
%   kb = element matrix for bending stiffness
%   ks = element matrix for shear stiffness
%   f = element vector
%   kk = system matrix                                             
%   ff = system vector                                                 
%   disp = system nodal displacement vector
%   gcoord = coordinate values of each node
%   nodes = nodal connectivity of each element
%   index = a vector containing system dofs associated with each element     
%   pointb = matrix containing sampling points for bending term
%   weightb = matrix containing weighting coefficients for bending term
%   points = matrix containing sampling points for shear term
%   weights = matrix containing weighting coefficients for shear term
%   bcdof = a vector containing dofs associated with boundary conditions     
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'                                              
%   kinmtpb = matrix for kinematic equation for bending
%   matmtpb = matrix for material property for bending
%   kinmtps = matrix for kinematic equation for shear
%   matmtps = matrix for material property for shear
%
%----------------------------------------------------------------------------            

%------------------------------------
%  input data for control parameters
%------------------------------------

clear
nel=4;                   % number of elements
nnel=4;                  % number of nodes per element
ndof=3;                  % number of dofs per node
nnode=9;                 % total number of nodes in system
sdof=nnode*ndof;         % total system dofs  
edof=nnel*ndof;          % degrees of freedom per element
emodule=30e6;            % elastic modulus
poisson=0.3;             % Poisson's ratio
t=0.1;                   % plate thickness
nglxb=2; nglyb=2;        % 2x2 Gauss-Legendre quadrature for bending 
nglb=nglxb*nglyb;        % number of sampling points per element for bending
nglxs=1; nglys=1;        % 1x1 Gauss-Legendre quadrature for shear 
ngls=nglxs*nglys;        % number of sampling points per element for shear

%---------------------------------------------
%  input data for nodal coordinate values
%  gcoord(i,j) where i->node no. and j->x or y
%---------------------------------------------

gcoord=[0.0  0.0; 2.5  0.0; 5.0  0.0; 
        0.0  2.5; 2.5  2.5; 5.0  2.5;
        0.0  5.0; 2.5  5.0; 5.0  5.0];

%---------------------------------------------------------
%  input data for nodal connectivity for each element
%  nodes(i,j) where i-> element no. and j-> connected nodes
%---------------------------------------------------------

nodes=[1 2 5 4; 2 3 6 5; 4 5 8 7; 5 6 9 8];

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------

bcdof=[1 2 3 4 6 7 9 11 12 16 20 21 23 25 26];  % constrained dofs
bcval=zeros(1,15);        % whose described values are 0 

%----------------------------------------------
%  initialization of matrices and vectors
%----------------------------------------------

ff=zeros(sdof,1);       % system force vector
kk=zeros(sdof,sdof);    % system matrix
disp=zeros(sdof,1);     % system displacement vector
index=zeros(edof,1);    % index vector
kinmtpb=zeros(3,edof);   % kinematic matrix for bending
matmtpb=zeros(3,3);      % constitutive matrix for bending
kinmtps=zeros(2,edof);   % kinematic matrix for shear
matmtps=zeros(2,2);      % constitutive matrix for shear

%----------------------------
%  force vector
%----------------------------

ff(27)=10;              % transverse force applied at node 9 

%-----------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%-----------------------------------------------------------------
%
%  for bending stiffness
%
[pointb,weightb]=feglqd2(nglxb,nglyb);     % sampling points & weights
matmtpb=fematiso(1,emodule,poisson)*t^3/12;  % bending material property
%
%  for shear stiffness
%
[points,weights]=feglqd2(nglxs,nglys);     % sampling points & weights
shearm=0.5*emodule/(1.0+poisson);          % shear modulus
shcof=5/6;                                 % shear correction factor
matmtps=shearm*shcof*t*[1 0; 0 1];         % shear material property

for iel=1:nel           % loop for the total number of elements

for i=1:nnel
nd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
xcoord(i)=gcoord(nd(i),1);  % extract x value of the node
ycoord(i)=gcoord(nd(i),2);  % extract y value of the node
end

k=zeros(edof,edof);         % initialization of element matrix to zero
kb=zeros(edof,edof);        % initialization of bending matrix to zero
ks=zeros(edof,edof);        % initialization of shear matrix to zero

%------------------------------------------------------
%  numerical integration for bending term
%------------------------------------------------------

for intx=1:nglxb
x=pointb(intx,1);                  % sampling point in x-axis
wtx=weightb(intx,1);               % weight in x-axis
for inty=1:nglyb
y=pointb(inty,2);                  % sampling point in y-axis
wty=weightb(inty,2) ;              % weight in y-axis

[shape,dhdr,dhds]=feisoq4(x,y);     % compute shape functions and
                                    % derivatives at sampling point

jacob2=fejacob2(nnel,dhdr,dhds,xcoord,ycoord);  % compute Jacobian

detjacob=det(jacob2);                 % determinant of Jacobian
invjacob=inv(jacob2);                 % inverse of Jacobian matrix

[dhdx,dhdy]=federiv2(nnel,dhdr,dhds,invjacob); % derivatives w.r.t.
                                               % physical coordinate

kinmtpb=fekinepb(nnel,dhdx,dhdy);          % bending kinematic matrix

%--------------------------------------------
%  compute bending element matrix
%--------------------------------------------

kb=kb+kinmtpb'*matmtpb*kinmtpb*wtx*wty*detjacob;

end
end                      % end of numerical integration loop for bending term

%------------------------------------------------------
%  numerical integration for bending term
%------------------------------------------------------

for intx=1:nglxs
x=points(intx,1);                  % sampling point in x-axis
wtx=weights(intx,1);               % weight in x-axis
for inty=1:nglys
y=points(inty,2);                  % sampling point in y-axis
wty=weights(inty,2) ;              % weight in y-axis

[shape,dhdr,dhds]=feisoq4(x,y);     % compute shape functions and
                                    % derivatives at sampling point

jacob2=fejacob2(nnel,dhdr,dhds,xcoord,ycoord);  % compute Jacobian

detjacob=det(jacob2);                 % determinant of Jacobian
invjacob=inv(jacob2);                 % inverse of Jacobian matrix

[dhdx,dhdy]=federiv2(nnel,dhdr,dhds,invjacob); % derivatives w.r.t.
                                               % physical coordinate

kinmtps=fekineps(nnel,dhdx,dhdy,shape);        % shear kinematic matrix

%----------------------------------------
%  compute shear element matrix
%----------------------------------------

ks=ks+kinmtps'*matmtps*kinmtps*wtx*wty*detjacob;   

end
end                      % end of numerical integration loop for shear term

%--------------------------------
%  compute element matrix
%--------------------------------

k=kb+ks;

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

