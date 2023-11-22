%----------------------------------------------------------------------------%
% Example 7.5.1                                                              
% to solve natural frequency of 1-d bar structure
%                                                                            
% Problem description                                                        
%   Find the natural frequency of a bar structure           
%   as shown in Fig. 7.5.1.           
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
nel=4;           % number of elements
nnel=2;          % number of nodes per element
ndof=1;          % number of dofs per node
nnode=5;         % total number of nodes in system
sdof=nnode*ndof; % total system dofs  

%---------------------------
%  nodal coordinates
%---------------------------

gcoord(1,1)=0.0;  
gcoord(2,1)=1.0; 
gcoord(3,1)=2.0;  
gcoord(4,1)=3.0;  
gcoord(5,1)=4.0;

%------------------------------------------
%  material and geometric properties
%------------------------------------------

prop(1)=200e9;       % elastic modulus 
prop(2)=0.001;       % cross-sectional area
prop(3)=7860;        % density

%-----------------------------
%  nodal connectivity
%-----------------------------

nodes(1,1)=1;  nodes(1,2)=2;   
nodes(2,1)=2;  nodes(2,2)=3;   
nodes(3,1)=3;  nodes(3,2)=4;   
nodes(4,1)=4;  nodes(4,2)=5;   

%-----------------------------
%  applied constraints
%-----------------------------

%bcdof(1)=1;      % 1st dof is constrained

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

x1=gcoord(nd(1),1);   % coordinate of 1st node
x2=gcoord(nd(2),1);   % coordinate of 2nd node

leng=(x2-x1);         % element length

el=prop(1);           % extract elastic modulus
area=prop(2);         % extract cross-sectional area
rho=prop(3);          % extract mass density

index=feeldof(nd,nnel,ndof);  % extract system dofs for the element

ipt=1;                        % flag for consistent mass matrix
[k,m]=fetruss1(el,leng,area,rho,ipt); % element matrix

kk=feasmbl1(kk,k,index);           % assemble system stiffness matrix
mm=feasmbl1(mm,m,index);           % assemble system mass matrix

end

%------------------------------
%  solve for eigenvalues
%------------------------------

%[kk,mm]=feaplycs(kk,mm,bcdof);  % apply the boundary conditions

fsol=eig(kk,mm);
fsol=sqrt(fsol);

%----------------------------
% print fem solutions
%----------------------------

num=1:1:sdof;
freqcy=[num' fsol]          % print natural frequency

%--------------------------------------------------------------------

