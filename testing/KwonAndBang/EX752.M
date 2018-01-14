%----------------------------------------------------------------------------%
% Example 7.5.2                                                              
% to solve natural frequency of 2-d truss structure
%                                                                            
% Problem description                                                        
%   Find the natural frequency of a truss structure           
%   as shown in Fig. 7.4.2.           
%                                                                            
% Variable descriptions                                                      
%   k = element stiffness matrix                                             
%   m = element mass matrix
%   kk = system stiffness matrix                                             
%   mm = system mass vector                                                 
%   index = a vector containing system dofs associated with each element     
%   gcoord = global coordinate matrix
%   prop = element property matrix
%   nodes = nodal connectivity matrix for each element
%   bcdof = a vector containing dofs associated with boundary conditions     
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'                                              
%----------------------------------------------------------------------------%            

%---------------------------
%  control input data
%---------------------------

clear
nel=9;           % number of elements
nnel=2;          % number of nodes per element
ndof=2;          % number of dofs per node
nnode=6;         % total number of nodes in system
sdof=nnode*ndof; % total system dofs  

%---------------------------
%  nodal coordinates
%---------------------------

gcoord(1,1)=0.0;  gcoord(1,2)=0.0;  
gcoord(2,1)=4.0;  gcoord(2,2)=0.0;   
gcoord(3,1)=4.0;  gcoord(3,2)=3.0;  
gcoord(4,1)=8.0;  gcoord(4,2)=0.0;  
gcoord(5,1)=8.0;  gcoord(5,2)=3.0;   
gcoord(6,1)=12.;  gcoord(6,2)=0.0;  

%------------------------------------------
%  material and geometric properties
%------------------------------------------

prop(1)=200e9;        % elastic modulus 
prop(2)=0.0025;       % cross-sectional area
prop(3)=7860;         % density

%-----------------------------
%  nodal connectivity
%-----------------------------

nodes(1,1)=1;  nodes(1,2)=2;   
nodes(2,1)=1;  nodes(2,2)=3;   
nodes(3,1)=2;  nodes(3,2)=3;   
nodes(4,1)=2;  nodes(4,2)=4;   
nodes(5,1)=3;  nodes(5,2)=4;   
nodes(6,1)=3;  nodes(6,2)=5;   
nodes(7,1)=4;  nodes(7,2)=5;   
nodes(8,1)=4;  nodes(8,2)=6;   
nodes(9,1)=5;  nodes(9,2)=6;   

%-----------------------------
%  applied constraints
%-----------------------------

bcdof(1)=1;      % 1st dof (horizontal displ) is constrained
bcval(1)=0;      % whose described value is 0 
bcdof(2)=2;      % 2nd dof (vertical displ) is constrained
bcdof(3)=12;     % 12th dof (horizontal displ) is constrained
bcval(3)=0;      % whose described value is 0 

%----------------------------
%  initialization to zero
%----------------------------

kk=zeros(sdof,sdof);           % system stiffness matrix
mm=zeros(sdof,sdof);           % system mass matrix
index=zeros(nnel*ndof,1);      % index vector

%--------------------------
%  loop for elements
%--------------------------

for iel=1:nel    % loop for the total number of elements

nd(1)=nodes(iel,1);   % 1st connected node for the (iel)-th element
nd(2)=nodes(iel,2);   % 2nd connected node for the (iel)-th element

x1=gcoord(nd(1),1); y1=gcoord(nd(1),2);  % coordinate of 1st node
x2=gcoord(nd(2),1); y2=gcoord(nd(2),2);  % coordinate of 2nd node

leng=sqrt((x2-x1)^2+(y2-y1)^2);  % element length

if (x2-x1)==0; 
if y2>y1;   
   beta=2*atan(1);       % angle between local and global axes
else
   beta=-2*atan(1);
end
else
beta=atan((y2-y1)/(x2-x1));
end

el=prop(1);           % extract elastic modulus
area=prop(2);         % extract cross-sectional area
rho=prop(3);          % extract mass density

index=feeldof(nd,nnel,ndof);  % extract system dofs for the element

ipt=1;                        % flag for consistent mass matrix
[k,m]=fetruss2(el,leng,area,rho,beta,ipt); % element matrix

kk=feasmbl1(kk,k,index);           % assemble system stiffness matrix
mm=feasmbl1(mm,m,index);           % assemble system mass matrix

end

%-------------------------------------------
%  apply constraints and solve 
%-------------------------------------------

[kk,mm]=feaplycs(kk,mm,bcdof);  % apply the boundary conditions

fsol=eig(kk,mm);
fsol=sqrt(fsol);

%----------------------------
% print fem solutions
%----------------------------

num=1:1:sdof;
freqcy=[num' fsol]          % print natural frequency

%--------------------------------------------------------------------

