%----------------------------------------------------------------------------
% Example 10.10.2                                                              
%   A cylinder is pinched across its diagonal direction with a force of 
%   100 lb at center lengthwise. It has the radius of 5 in., length of
%   10.35 in., and thickness of 0.094 in. 
%   The material has elastic modulus of 10.5 msi and Poisson's ratio 0.3125.
%   A one-eigths of the cylinder is modeled due to symmetry using          
%   25 four-node shell elements. 
%   (see Fig. 10.10.2 for the finite element mesh)
%
% Variable descriptions                                                      
%   k = element matrix in terms of local axes
%   kt = element matrix in terms of global axes
%   f = element vector
%   kk = system matrix                                             
%   ff = system vector                                                 
%   disp = system nodal displacement vector
%   gcoord = coordinate values of each node
%   nodes = nodal connectivity of each element
%   index = a vector containing system dofs associated with each element     
%   point1 = matrix containing sampling points for inplne axes
%   weight1 = matrix containing weighting coefficients for inplane axes
%   pointz = matrix containing sampling points for transverse axis
%   weightz = matrix containing weighting coefficients for transverse axis
%   bcdof = a vector containing dofs associated with boundary conditions     
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'                                              
%   bmtx = matrix for kinematic equation 
%   dmtx = matrix for material property 
%   trsh = nodal variable transformation matrix
%   rot6 = strain transformation matrix
%   aj = 3-D Jacobian matrix
%
%----------------------------------------------------------------------------            

%------------------------------------
%  input data for control parameters
%------------------------------------

clear
nel=25;                   % number of elements
nnel=4;                  % number of nodes per element
ndof=6;                  % number of dofs per node
nnode=36;                 % total number of nodes in system
sdof=nnode*ndof;         % total system dofs  
edof=nnel*ndof;          % degrees of freedom per element
emodule=10.5e6;            % elastic modulus
poisson=0.3125;             % Poisson's ratio
thick=0.094;                   % shell thickness
nglx=1; ngly=1; nglz=2;  % 1x1x2 Gauss-Legendre quadrature  
ngl=nglx*ngly*nglz;   % number of sampling points per element 

%---------------------------------------------
%  input data for nodal coordinate values
%  gcoord(i,j) where i->node no. and j->x,y,z
%---------------------------------------------

gcoord=[
    0.0000    0.0000    5.0000
    0.0000    1.0350    5.0000
    0.0000    2.0700    5.0000
    0.0000    3.1050    5.0000
    0.0000    4.1400    5.0000
    0.0000    5.1750    5.0000
    1.5451    0.0000    4.7553
    1.5451    1.0350    4.7553
    1.5451    2.0700    4.7553
    1.5451    3.1050    4.7553
    1.5451    4.1400    4.7553
    1.5451    5.1750    4.7553
    2.9389    0.0000    4.0451
    2.9389    1.0350    4.0451
    2.9389    2.0700    4.0451
    2.9389    3.1050    4.0451
    2.9389    4.1400    4.0451
    2.9389    5.1750    4.0451
    4.0451    0.0000    2.9389
    4.0451    1.0350    2.9389
    4.0451    2.0700    2.9389
    4.0451    3.1050    2.9389
    4.0451    4.1400    2.9389
    4.0451    5.1750    2.9389
    4.7553    0.0000    1.5451
    4.7553    1.0350    1.5451
    4.7553    2.0700    1.5451
    4.7553    3.1050    1.5451
    4.7553    4.1400    1.5451
    4.7553    5.1750    1.5451
    5.0000    0.0000    0.0000
    5.0000    1.0350    0.0000
    5.0000    2.0700    0.0000
    5.0000    3.1050    0.0000
    5.0000    4.1400    0.0000
    5.0000    5.1750    0.0000
];
     
%---------------------------------------------------------
%  input data for nodal connectivity for each element
%  nodes(i,j) where i-> element no. and j-> connected nodes
%---------------------------------------------------------

nodes=[
 7    8    2    1  ;
 8    9    3    2  ;
 9   10    4    3  ;
10   11    5    4  ;
11   12    6    5  ;
13   14    8    7  ;
14   15    9    8  ;
15   16   10    9  ;
16   17   11   10  ;
17   18   12   11  ;
19   20   14   13  ;
20   21   15   14  ;
21   22   16   15  ;
22   23   17   16  ;
23   24   18   17  ;
25   26   20   19  ;
26   27   21   20  ;
27   28   22   21  ;
28   29   23   22  ;
29   30   24   23  ;
31   32   26   25  ;
32   33   27   26  ;
33   34   28   27  ;
34   35   29   28  ;
35   36   30   29  ];

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------

bcdof=[1               5   6  ...
       7              11  12  ...
      13              17  18  ...
      19              23  24  ... 
      25              29  30  ...
      31  32      34  35  36   ...
          38                   ...
          68      70      72   ...
         104     106     108  ...
         140     142     144  ...
         176     178     180  ...
             183 184 185      ...
             189 190 191      ...
             195 196 197      ...
             201 202 203      ...
             207 208 209      ...
         212 213 214 215 216 ];        % constrained dofs
bcval=zeros(size(bcdof));        % whose described values are 0 

%----------------------------------------------
%  initialization of matrices and vectors
%----------------------------------------------

ff=zeros(sdof,1);       % system force vector
kk=zeros(sdof,sdof);    % system matrix
disp=zeros(sdof,1);     % system displacement vector
index=zeros(edof,1);    % index vector
dmtx=zeros(6,6);      % constitutive matrix for shear

%----------------------------
%  force vector
%----------------------------

ff(33)=-25;              % transverse force at node 3

%-----------------------------------------------------------------
%  compute material property matrix
%-----------------------------------------------------------------

dmtx(1,1)=emodule/(1.0-poisson^2);
dmtx(1,2)=poisson*dmtx(1,1);
dmtx(2,1)=dmtx(1,2);
dmtx(2,2)=dmtx(1,1);
dmtx(4,4)=emodule/(2.0*(1.0+poisson));
dmtx(5,5)=(5/6)*dmtx(4,4);
dmtx(6,6)=dmtx(5,5);

%-----------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%-----------------------------------------------------------------
%
%  for inplane integration
%
[point1,weight1]=feglqd2(nglx,ngly);     % sampling points & weights
%
%  for transverse integration
%
[pointz,weightz]=feglqd1(nglz);     % sampling points & weights

for iel=1:nel           % loop for the total number of elements

for i=1:nnel
nd(i)=nodes(iel,i);         % extract connected node for (iel)-th element
xcoord(i)=gcoord(nd(i),1);  % extract x value of the node
ycoord(i)=gcoord(nd(i),2);  % extract y value of the node
zcoord(i)=gcoord(nd(i),3);  % extract z value of the node
end

k=zeros(edof,edof);         % initialization of element matrix
kt=zeros(edof,edof);        % initialization of transformed element matrix
trsh=zeros(edof,edof);      % initialization of coordinate transformation matrix

%------------------------------------------------------
%  compute direction cosine vectors
%  v1 and v2: tangent to the element
%  v3: normal to the element
%------------------------------------------------------

d1(1)=xcoord(2)-xcoord(1);		% define the vector from node 1 to node 2
d1(2)=ycoord(2)-ycoord(1);
d1(3)=zcoord(2)-zcoord(1);

d2(1)=xcoord(4)-xcoord(1);		% define the vector from node 4 to node 1
d2(2)=ycoord(4)-ycoord(1);
d2(3)=zcoord(4)-zcoord(1);

v3(1)=d1(2)*d2(3)-d1(3)*d2(2); % vector v3 normal to vectors d1 and d2
v3(2)=d1(3)*d2(1)-d1(1)*d2(3);
v3(3)=d1(1)*d2(2)-d1(2)*d2(1);
sum=sqrt(v3(1)^2+v3(2)^2+v3(3)^2); 
v3(1)=v3(1)/sum;	% make vector v3 a unit vector
v3(2)=v3(2)/sum;
v3(3)=v3(3)/sum;

if v3(2) > 0.999999	 % if v3 is along the y-axis, v1 is set along x-axis
v1(1)=1.0;
v1(2)=0.0;
v1(3)=0.0;
else					% otherwise, v1 is the cross product of y-axis and v3
v1(1)=v3(3);
v1(2)=0.0;
v1(3)=-v3(1);
sum=sqrt(v1(1)^2+v1(2)^2+v1(3)^2); 
v1(1)=v1(1)/sum;	% make vector v1 a unit vector
v1(2)=v1(2)/sum;
v1(3)=v1(3)/sum;
end

v2(1)=v3(2)*v1(3)-v3(3)*v1(2);   % v2 is cross product of v3 and v1
v2(2)=v3(3)*v1(1)-v3(1)*v1(3);
v2(3)=v3(1)*v1(2)-v3(2)*v1(1);
sum=sqrt(v2(1)^2+v2(2)^2+v2(3)^2); 
v2(1)=v2(1)/sum;	% make vector v2 a unit vector
v2(2)=v2(2)/sum;
v2(3)=v2(3)/sum;

%---------------------------------------------------------
%   construct nodal variable transformation matrix 
%---------------------------------------------------------

a(1,1)=v1(1);
a(2,1)=v1(2);
a(3,1)=v1(3);
a(1,2)=v2(1);
a(2,2)=v2(2);
a(3,2)=v2(3);
a(1,3)=v3(1);
a(2,3)=v3(2);
a(3,3)=v3(3);

ainv=inv(a);

for i=1:nnel
i1=(i-1)*ndof+1;
i2=i1+1;
i3=i2+1;
i4=i3+1;
i5=i4+1;
i6=i5+1;

trsh(i1,i1)=1.0;
trsh(i2,i2)=1.0;
trsh(i3,i3)=1.0;
trsh(i4,i4)=ainv(1,1);
trsh(i4,i5)=ainv(1,2);
trsh(i4,i6)=ainv(1,3);
trsh(i5,i4)=ainv(2,1);
trsh(i5,i5)=ainv(2,2);
trsh(i5,i6)=ainv(2,3);
trsh(i6,i4)=ainv(3,1);
trsh(i6,i5)=ainv(3,2);
trsh(i6,i6)=ainv(3,3);
end
  
%------------------------------------------------------
%  strain transformation matrix
%------------------------------------------------------

rot6(1,1)=v1(1)^2;
rot6(1,2)=v1(2)^2;
rot6(1,3)=v1(3)^2;
rot6(1,4)=v1(1)*v1(2);
rot6(1,5)=v1(2)*v1(3);
rot6(1,6)=v1(1)*v1(3);
rot6(2,1)=v2(1)^2;
rot6(2,2)=v2(2)^2;
rot6(2,3)=v2(3)^2;
rot6(2,4)=v2(1)*v2(2);
rot6(2,5)=v2(2)*v2(3);
rot6(2,6)=v2(1)*v2(3);
rot6(3,1)=v3(1)^2;
rot6(3,2)=v3(2)^2;
rot6(3,3)=v3(3)^2;
rot6(3,4)=v3(1)*v3(2);
rot6(3,5)=v3(2)*v3(3);
rot6(3,6)=v3(1)*v3(3);
rot6(4,1)=2.0*v1(1)*v2(1);
rot6(4,2)=2.0*v1(2)*v2(2);
rot6(4,3)=2.0*v1(3)*v2(3);
rot6(4,4)=v1(1)*v2(2)+v2(1)*v1(2);
rot6(4,5)=v1(2)*v2(3)+v2(2)*v1(3);
rot6(4,6)=v1(3)*v2(1)+v2(3)*v1(1);
rot6(5,1)=2.0*v2(1)*v3(1);
rot6(5,2)=2.0*v2(2)*v3(2);
rot6(5,3)=2.0*v2(3)*v3(3);
rot6(5,4)=v2(1)*v3(2)+v3(1)*v2(2);
rot6(5,5)=v2(2)*v3(3)+v3(2)*v2(3);
rot6(5,6)=v2(3)*v3(1)+v3(3)*v2(1);
rot6(6,1)=2.0*v3(1)*v1(1);
rot6(6,2)=2.0*v3(2)*v1(2);
rot6(6,3)=2.0*v3(3)*v1(3);
rot6(6,4)=v3(1)*v1(2)+v1(1)*v3(2);
rot6(6,5)=v3(2)*v1(3)+v1(2)*v3(3);
rot6(6,6)=v3(3)*v1(1)+v1(3)*v3(1);

%------------------------------------------------------
%  material property matrix transformation
%------------------------------------------------------

dmtxt=rot6'*dmtx*rot6;

%------------------------------------------------------
%  numerical integration loop
%------------------------------------------------------

for intx=1:nglx
r=point1(intx);                    % sampling point in x-axis
wtx=weight1(intx,1);               % weight in x-axis
for inty=1:ngly
s=point1(inty,2);                  % sampling point in y-axis
wty=weight1(inty,2) ;              % weight in y-axis
for intz=1:nglz
t=pointz(intz);                  % sampling point in z-axis
wtz=weightz(intz) ;              % weight in z-axis

[shape2,dhdr,dhds]=feisoq4(r,s);     % compute 2-D shape functions and
                                    % derivatives at sampling point

%[shape1,dhdz]=feisol2(t);     % compute 1-D shape functions and
                                    % derivatives at sampling point
shape1(1)=0.5*(1-t);
shape1(2)=0.5*(1+t);
dhdt(1)=-0.5;
dhdt(2)=0.5;

%--------------------------------------------------
%   compute jacobian matrix and its inverse
%--------------------------------------------------

hz=thick*0.5*t;
hzdt=thick*0.5;

aj(1,1)=dhdr(1)*xcoord(1)+dhdr(2)*xcoord(2)+dhdr(3)*xcoord(3)+ ...
   dhdr(4)*xcoord(4)+dhdr(1)*hz*v3(1)+dhdr(2)*hz*v3(1)+ ...
   dhdr(3)*hz*v3(1)+dhdr(4)*hz*v3(1);
aj(2,1)=dhds(1)*xcoord(1)+dhds(2)*xcoord(2)+dhds(3)*xcoord(3)+ ...
   dhds(4)*xcoord(4)+dhds(1)*hz*v3(1)+dhds(2)*hz*v3(1)+ ...
   dhds(3)*hz*v3(1)+dhds(4)*hz*v3(1);
aj(3,1)=hzdt*v3(1);
aj(1,2)=dhdr(1)*ycoord(1)+dhdr(2)*ycoord(2)+dhdr(3)*ycoord(3)+ ...
   dhdr(4)*ycoord(4)+dhdr(1)*hz*v3(2)+dhdr(2)*hz*v3(2)+ ...
   dhdr(3)*hz*v3(2)+dhdr(4)*hz*v3(2);
aj(2,2)=dhds(1)*ycoord(1)+dhds(2)*ycoord(2)+dhds(3)*ycoord(3)+ ...
   dhds(4)*ycoord(4)+dhds(1)*hz*v3(2)+dhds(2)*hz*v3(2)+ ...
   dhds(3)*hz*v3(2)+dhds(4)*hz*v3(2);
aj(3,2)=hzdt*v3(2);
aj(1,3)=dhdr(1)*zcoord(1)+dhdr(2)*zcoord(2)+dhdr(3)*zcoord(3)+ ...
   dhdr(4)*zcoord(4)+dhdr(1)*hz*v3(3)+dhdr(2)*hz*v3(3)+ ...
   dhdr(3)*hz*v3(3)+dhdr(4)*hz*v3(3);
aj(2,3)=dhds(1)*zcoord(1)+dhds(2)*zcoord(2)+dhds(3)*zcoord(3)+ ...
   dhds(4)*zcoord(4)+dhds(1)*hz*v3(3)+dhds(2)*hz*v3(3)+ ...
   dhds(3)*hz*v3(3)+dhds(4)*hz*v3(3);
aj(3,3)=hzdt*v3(3);

ajinv=inv(aj);
det3=det(aj);

%-----------------------------------------------------------
%   compute global derivatives
%-----------------------------------------------------------

for i=1:nnel
derivg(1,i)=ajinv(1,1)*dhdr(i)+ajinv(1,2)*dhds(i);
derivg(2,i)=ajinv(2,1)*dhdr(i)+ajinv(2,2)*dhds(i);
derivg(3,i)=ajinv(3,1)*dhdr(i)+ajinv(3,2)*dhds(i);
dhdz(1,i)=ajinv(1,3)*hzdt;
dhdz(2,i)=ajinv(2,3)*hzdt;
dhdz(3,i)=ajinv(3,3)*hzdt;
end

%-----------------------------------------------------------
%   compute strain-displacement matrix called bmtx
%-----------------------------------------------------------

bmtx=zeros(ndof,edof);

for i=1:nnel
i1=(i-1)*ndof+1;
i2=i1+1;
i3=i2+1;
i4=i3+1;
i5=i4+1;
i6=i5+1;

gk1=derivg(1,i)*hz+shape2(i)*dhdz(1,i);
gk2=derivg(2,i)*hz+shape2(i)*dhdz(2,i);
gk3=derivg(3,i)*hz+shape2(i)*dhdz(3,i);

bmtx(1,i1)=derivg(1,i);	% elements associated with epsilon_x
bmtx(1,i4)=gk1*(-v2(1)); 
bmtx(1,i5)=gk1*v1(1);
bmtx(2,i2)=derivg(2,i); % elements related to epsilon_y
bmtx(2,i4)=gk2*(-v2(2));
bmtx(2,i5)=gk2*v1(2);
bmtx(3,i3)=derivg(3,i); % elements related to epsilon_z
bmtx(3,i4)=gk3*(-v2(3));
bmtx(3,i5)=gk3*v1(3);
bmtx(4,i1)=derivg(2,i); % elements related to gamma_xy
bmtx(4,i2)=derivg(1,i);
bmtx(4,i4)=gk2*(-v2(1))+gk1*(-v2(2)); 
bmtx(4,i5)=gk2*v1(1)+gk1*v1(2);
bmtx(5,i2)=derivg(3,i); % elements related to gamma_yz
bmtx(5,i3)=derivg(2,i);
bmtx(5,i4)=gk3*(-v2(2))+gk2*(-v2(3)); 
bmtx(5,i5)=gk3*v1(2)+gk2*v1(3);
bmtx(6,i1)=derivg(3,i); % elements related to gamma_xz
bmtx(6,i3)=derivg(1,i);
bmtx(6,i4)=gk3*(-v2(1))+gk1*(-v2(3)); 
bmtx(6,i5)=gk3*v1(1)+gk1*v1(3);
end

detwt=det3*wtx*wty*wtz;

k=k+bmtx'*dmtxt*bmtx*detwt;  % element matrix in local axes

end
end
end	% end of numerical integration loop
   
kt=trsh'*k*trsh;   % element matrix in global axis


index=feeldof(nd,nnel,ndof);% extract system dofs associated with element

kk=feasmbl1(kk,kt,index);  % assemble element matrices 

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
displace=[num' disp]                % print nodal displacements
     			     
%------------------------------------------------------------------------

