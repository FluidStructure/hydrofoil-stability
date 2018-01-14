%----------------------------------------------------------------------------%
% Example 8.10.3                                                             %
% to find the natural frequencies for a 2-d frame using frame elements       %
%                                                                            %
% Problem description                                                        %
%   Find the natural frequencies of a frame of L-shape which is made of      %
%   two beams of length of 1 m each. Both beams have                         %
%   cross-sections of 0.01 m by 0.01 m. The elastic modulus is 100 GPa.      %
%   The beam has mass density of 1000 Kg/m^3.  Use 10 elements.              %
%                                                                            %
% Variable descriptions                                                      %
%   x and y = global x and y coordiates of each node                         %
%   k = element stiffness matrix                                             %
%   kk = system stiffness matrix                                             %
%   m = element mass matrix                                                  %
%   mm = system mass matrix                                                  %
%   index = a vector containing system dofs associated with each element     %
%   bcdof = a vector containing dofs associated with boundary conditions     %
%----------------------------------------------------------------------------%            

clear
nel=10;          % number of elements
nnel=2;          % number of nodes per element
ndof=3;          % number of dofs per node
nnode=(nnel-1)*nel+1;   % total number of nodes in system
sdof=nnode*ndof; % total system dofs  

x(1)=0; y(1)=0;     % x, y coord. values of node 1 in terms of the global axis
x(2)=0; y(2)=0.2;   % x, y coord. values of node 2 in terms of the global axis  
x(3)=0; y(3)=0.4;   % x, y coord. values of node 3 in terms of the global axis
x(4)=0; y(4)=0.6;   % x, y coord. values of node 4 in terms of the global axis  
x(5)=0; y(5)=0.8;   % x, y coord. values of node 5 in terms of the global axis
x(6)=0; y(6)=1;     % x, y coord. values of node 6 in terms of the global axis  
x(7)=0.2; y(7)=1;   % x, y coord. values of node 7 in terms of the global axis
x(8)=0.4; y(8)=1;   % x, y coord. values of node 8 in terms of the global axis
x(9)=0.6; y(9)=1;   % x, y coord. values of node 9 in terms of the global axis  
x(10)=0.8; y(10)=1; % x, y coord. values of node 10 in terms of the global axis
x(11)=1; y(11)=1;   % x, y coord. values of node 11 in terms of the global axis  

el=100*10^9;        % elastic modulus
area=0.0001;        % cross-sectional area
xi=8.3333*10^(-10); % moment of inertia of cross-section
rho=1000;           % mass density per volume (dummy value for static analysis)
bcdof(1)=1;         % transverse deflection at node 1 is constrained
bcdof(2)=2;         % axial displacement at node 1 is constrained
bcdof(3)=3;         % slope at node 1 is constrained

kk=zeros(sdof,sdof);      % initialization of system stiffness matrix
mm=zeros(sdof,sdof);      % initialization of system mass matrix
index=zeros(nel*ndof,1);  % initialization of index vector

for iel=1:nel    % loop for the total number of elements

index=feeldof1(iel,nnel,ndof);  % extract system dofs associated with element

node1=iel;      % starting node number for element 'iel'
node2=iel+1;    % ending node number for element 'iel'

x1=x(node1); y1=y(node1); % x and y coordinate values of 'node1'
x2=x(node2); y2=y(node2); % x and y coordinate values of 'node2'

leng=sqrt((x2-x1)^2+(y2-y1)^2); % length of element 'iel'

if (x2-x1)==0;  % compute the angle between the local and global axes    
   beta=pi/2; 
else
   beta=atan((y2-y1)/(x2-x1));
end 

[k,m]=feframe2(el,xi,leng,area,rho,beta,1); % compute element stiffness matrix

kk=feasmbl1(kk,k,index); % assemble element matrices into system matrix

mm=feasmbl1(mm,m,index); % assemble element mass matrices into system matrix

end

[kn,mn]=feaplycs(kk,mm,bcdof);  % apply the boundary conditions

fsol=eig(kn,mn);   % solve the matrix equation and print
fsol=sqrt(fsol)


%_______________________________________________________________________
