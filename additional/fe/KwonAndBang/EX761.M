%----------------------------------------------------------------------------%
% Example 7.6.1                                                              
% to solve transient response of 1-d bar structure
%                                                                            
% Problem description                                                        
%   Find the dynamic behavior of a bar structure,           
%   as shown in Fig. 7.6.1, subjected to a step 
%   force function at the right end.           
%                                                                            
% Variable descriptions                                                      
%   k = element stiffness matrix                                             
%   m = element mass matrix
%   kk = system stiffness matrix                                             
%   mm = system mass vector                                                 
%   ff = system force vector
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
nel=10;             % number of elements
nnel=2;             % number of nodes per element
ndof=1;             % number of dofs per node
nnode=11;           % total number of nodes in system
sdof=nnode*ndof;    % total system dofs  
dt=0.0001;          % time step size
ti=0;               % initial time
tf=0.05;           % final time
nt=fix((tf-ti)/dt); % number of time steps

%---------------------------
%  nodal coordinates
%---------------------------

gcoord(1,1)=0.0;  
gcoord(2,1)=1.0; 
gcoord(3,1)=2.0;  
gcoord(4,1)=3.0;  
gcoord(5,1)=4.0;
gcoord(6,1)=5.0;  
gcoord(7,1)=6.0; 
gcoord(8,1)=7.0;  
gcoord(9,1)=8.0;  
gcoord(10,1)=9.0;
gcoord(11,1)=10.0;

%------------------------------------------
%  material and geometric properties
%------------------------------------------

prop(1)=200e9;       % elastic modulus 
prop(2)=0.001;       % cross-sectional area
prop(3)=7860;        % density

%-----------------------------
%  nodal connectivity
%-----------------------------

nodes(1,1)=1;   nodes(1,2)=2;   
nodes(2,1)=2;   nodes(2,2)=3;   
nodes(3,1)=3;   nodes(3,2)=4;   
nodes(4,1)=4;   nodes(4,2)=5;   
nodes(5,1)=5;   nodes(5,2)=6;   
nodes(6,1)=6;   nodes(6,2)=7;   
nodes(7,1)=7;   nodes(7,2)=8;   
nodes(8,1)=8;   nodes(8,2)=9;   
nodes(9,1)=9;   nodes(9,2)=10;   
nodes(10,1)=10; nodes(10,2)=11;   

%-----------------------------
%  applied constraints
%-----------------------------

nbc=1;           % number of constraints
bcdof(1)=1;      % 1st dof is constrained

%----------------------------
%  initialization to zero
%----------------------------

kk=zeros(sdof,sdof);           % system stiffness matrix
mm=zeros(sdof,sdof);           % system mass matrix
ff=zeros(sdof,1);              % system force vector
index=zeros(nnel*ndof,1);      % index vector
acc=zeros(sdof,nt);            % acceleartion matrix
vel=zeros(sdof,nt);            % velocity matrix
disp=zeros(sdof,nt);           % displacement matrix

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

%-------------------------------
%  initial condition
%-------------------------------

vel(:,1)=zeros(sdof,1);         % initial zero velocity
disp(:,1)=zeros(sdof,1);        % initial zero displacement

ff(11)=200;                     % step force at node 11 

%--------------------------------------------------------
% central difference scheme for time integration
%--------------------------------------------------------                                                        

mm=inv(mm);        % invert the mass matrix

for it=1:nt

acc(:,it)=mm*(ff-kk*disp(:,it));

  for i=1:nbc
  ibc=bcdof(i);  
  acc(ibc,it)=0;
  end

vel(:,it+1)=vel(:,it)+acc(:,it)*dt;
disp(:,it+1)=disp(:,it)+vel(:,it+1)*dt;

end

acc(:,nt+1)=mm*(ff-kk*disp(:,nt+1));

time=0:dt:nt*dt;
plot(time,disp(11,:))
xlabel('Time(seconds)')
ylabel('Tip displ. (m)')

%---------------------------------------------------------------------------









































