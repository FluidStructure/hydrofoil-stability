%----------------------------------------------------------------------------
% Example
%   A barrel vault has radius r=25 ft, subtended angle of 80 degrees,
%   length 50 f, and thickness 3in. 
%   Elastic modulus is 3 msi, Poisson's ratio is 0, and
%   weight of the shell is 90 lb per square ft. Two curved edges are
%   supported by rigid diaphragm and the other two edges are free. 
%   Solve using 4 by 4 elements of a quarter of the shell.         
%
% Variable descriptions                                                      
%   k = element matrix in the local axes                                            
%   ke = element matrix in the global axes
%   kb = element matrix for bending stiffness
%   ks = element matrix for shear stiffness
%   km = element matrix for membrane stiffness
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
%   kinmtsb = matrix for kinematic equation for bending
%   matmtsb = matrix for material property for bending
%   kinmtss = matrix for kinematic equation for shear
%   matmtss = matrix for material property for shear
%   kinmtsm = matrix for kinematic equation for membrane
%   matmtsm = matrix for material property for membrane
%   tr3d = transformation matrix from local to global axes
%
%----------------------------------------------------------------------------            

%------------------------------------
%  input data for control parameters
%------------------------------------

clear
nel=16;                   % number of elements
nnel=4;                  % number of nodes per element
ndof=6;                  % number of dofs per node
nnode=25;                 % total number of nodes in system
sdof=nnode*ndof;         % total system dofs  
edof=nnel*ndof;          % degrees of freedom per element
emodule=3e6;            % elastic modulus
poisson=0.0;             % Poisson's ratio
t=3;                   % plate thickness
nglxb=2; nglyb=2;        % 2x2 Gauss-Legendre quadrature for bending 
nglb=nglxb*nglyb;        % number of sampling points per element for bending
nglxs=1; nglys=1;        % 1x1 Gauss-Legendre quadrature for shear 
ngls=nglxs*nglys;        % number of sampling points per element for shear

%---------------------------------------------
%  input data for nodal coordinate values
%  gcoord(i,j) where i->node no. and j->x or y
%---------------------------------------------

gcoord=[0.0 0.0 0.0; 0.0 6.25 0.0; 0.0 12.5 0.0; 0.0 18.75 0.0;0.0 25.0 0.0; 
4.34 0.0 -0.38;4.34 6.25 -0.38;4.34 12.5 -0.38;4.34 18.75 -0.38;4.34 25.0 -0.38;     
8.55 0.0 -1.51;8.55 6.25 -1.51;8.55 12.5 -1.51;8.55 18.75 -1.51;8.55 25.0 -1.51;     
12.5 0.0 -3.35;12.5 6.25 -3.35;12.5 12.5 -3.35;12.5 18.75 -3.35;12.5 25.0 -3.35;     
16.1 0.0 -5.85;16.1 6.25 -5.85;16.1 12.5 -5.85;16.1 18.75 -5.85;16.1 25.0 -5.85];    
gcoord=12.0*gcoord;

%---------------------------------------------------------
%  input data for nodal connectivity for each element
%  nodes(i,j) where i-> element no. and j-> connected nodes
%---------------------------------------------------------

nodes=[6 7 2 1; 7 8 3 2; 8 9 4 3; 9 10 5 4;
      11 12 7 6; 12 13 8 7; 13 14 9 8; 14 15 10 9;
      16 17 12 11; 17 18 13 12; 18 19 14 13; 19 20 15 14;
      21 22 17 16; 22 23 18 17; 23 24 19 18; 24 25 20 19];

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------

bcdof=[1 3 5 6   7 11 12  13 17 18   19 23 24   25 26 28 29 30 ...
       31 33   56 58 60   61 63   86 88 90   91 93   116 118 120 ...
       121 122 123   146 148 150]; % constrained dofs
bcval=zeros(size(bcdof));     % whose described values are 0 

%----------------------------------------------
%  initialization of matrices and vectors
%----------------------------------------------

ff=zeros(sdof,1);        % system force vector
kk=zeros(sdof,sdof);     % system matrix
disp=zeros(sdof,1);      % system displacement vector
index=zeros(edof,1);     % index vector
kinmtsb=zeros(3,edof);   % kinematic matrix for bending
matmtsb=zeros(3,3);      % constitutive matrix for bending
kinmtsm=zeros(3,edof);   % kinematic matrix for membrane
matmtsm=zeros(3,3);      % constitutive matrix for membrane
kinmtss=zeros(2,edof);   % kinematic matrix for shear
matmtss=zeros(2,2);      % constitutive matrix for shear
tr3d=zeros(edof,edof);   % transformation matrix

%----------------------------
%  force vector
%----------------------------

ff(3)=-613.6;              % transverse force at node 1
ff(9)=-1227.2;              % transverse force at node 2
ff(15)=-1227.2;             % transverse force at node 3
ff(21)=-1227.2;              % transverse force at node 4
ff(27)=-613.6;              % transverse force at node 5
ff(33)=-1227.2;             % transverse force at node 6
ff(39)=-2454.4;              % transverse force at node 7
ff(45)=-2454.4;              % transverse force at node 8
ff(51)=-2454.4;             % transverse force at node 9
ff(57)=-1227.2;              % transverse force at node 10
ff(63)=-1227.2;              % transverse force at node 11
ff(69)=-2454.4;             % transverse force at node 12
ff(75)=-2454.4;              % transverse force at node 13
ff(81)=-2454.4;              % transverse force at node 14
ff(87)=-1227.2;             % transverse force at node 15
ff(93)=-1227.2;              % transverse force at node 16
ff(99)=-2454.4;              % transverse force at node 17
ff(105)=-2454.4;             % transverse force at node 18
ff(111)=-2454.4;              % transverse force at node 19
ff(117)=-1227.2;             % transverse force at node 20
ff(123)=-613.6;              % transverse force at node 21
ff(129)=-1227.2;              % transverse force at node 22
ff(135)=-1227.2;             % transverse force at node 23
ff(141)=-1227.2;              % transverse force at node 24
ff(147)=-613.6;             % transverse force at node 25

%-----------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%-----------------------------------------------------------------
%
%  for bending and membrane stiffness
%
[pointb,weightb]=feglqd2(nglxb,nglyb);     % sampling points & weights
matmtsm=fematiso(1,emodule,poisson)*t;       % membrane material property
matmtsb=fematiso(1,emodule,poisson)*t^3/12;  % bending material property
%%
%  for shear stiffness
%
[points,weights]=feglqd2(nglxs,nglys);     % sampling points & weights
shearm=0.5*emodule/(1.0+poisson);          % shear modulus
shcof=5/6;                                 % shear correction factor
matmtss=shearm*shcof*t*[1 0; 0 1];         % shear material property

for iel=1:nel           % loop for the total number of elements

for i=1:nnel
nd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
xcoord(i)=gcoord(nd(i),1);  % extract x value of the node
ycoord(i)=gcoord(nd(i),2);  % extract y value of the node
zcoord(i)=gcoord(nd(i),3);  % extract z value of the node
end

% 
%  compute the local direction cosines and local axes
%
[tr3d,xprime,yprime]=fetransh(xcoord,ycoord,zcoord,nnel);
%

k=zeros(edof,edof);         % element matrix in local axes
ke=zeros(edof,edof);        % element matrix in global axes
km=zeros(edof,edof);        % element membrane matrix
kb=zeros(edof,edof);        % element bending matrix
ks=zeros(edof,edof);        % element shear matrix

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

jacob2=fejacob2(nnel,dhdr,dhds,xprime,yprime);  % compute Jacobian

detjacob=det(jacob2);                 % determinant of Jacobian
invjacob=inv(jacob2);                 % inverse of Jacobian matrix

[dhdx,dhdy]=federiv2(nnel,dhdr,dhds,invjacob); % derivatives w.r.t.
                                               % physical coordinate

kinmtsb=fekinesb(nnel,dhdx,dhdy);          % bending kinematic matrix
kinmtsm=fekinesm(nnel,dhdx,dhdy);          % membrane kinematic matrix

%--------------------------------------------
%  compute bending element matrix
%--------------------------------------------

kb=kb+kinmtsb'*matmtsb*kinmtsb*wtx*wty*detjacob;
km=km+kinmtsm'*matmtsm*kinmtsm*wtx*wty*detjacob;
end
end                   % end of numerical integration loop for bending term

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

jacob2=fejacob2(nnel,dhdr,dhds,xprime,yprime);  % compute Jacobian

detjacob=det(jacob2);                 % determinant of Jacobian
invjacob=inv(jacob2);                 % inverse of Jacobian matrix
   
[dhdx,dhdy]=federiv2(nnel,dhdr,dhds,invjacob); % derivatives w.r.t.
                                               % physical coordinate

kinmtss=fekiness(nnel,dhdx,dhdy,shape);        % shear kinematic matrix

%----------------------------------------
%  compute shear element matrix
%----------------------------------------

ks=ks+kinmtss'*matmtss*kinmtss*wtx*wty*detjacob;   

end
end                      % end of numerical integration loop for shear term

%--------------------------------
%  compute element matrix
%--------------------------------

k=km+kb+ks;

%-----------------------------------------------
%  transform from local to global systems
%-----------------------------------------------

ke=tr3d'*k*tr3d;

index=feeldof(nd,nnel,ndof);% extract system dofs associated with element

kk=feasmbl1(kk,ke,index);  % assemble element matrices 

end

%-------------------------------------
%   check the singular drilling dof
%-------------------------------------

for i=1:sdof
if(abs(kk(i,i)) < 1e-5) 

sum=0.0;
for j=1:sdof
sum=sum+abs(kk(i,j));
end

if (sum < 1e-5) 
kk(i,i)=1;
end

end

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

%--------------------------------------------------------------

