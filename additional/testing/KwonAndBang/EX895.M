%----------------------------------------------------------------------------%
% Example 8.9.5                                                              %
% to solve a static beam deflection for a 2-d frame using frame elements     %
%                                                                            %
% Problem description                                                        %
%   Find the deflection of a frame of L-shape which is made of two beams     %
%   of lengths of 60 in. and 20 in., respectively. Both beams have           %
%   cross-sections of 2 in. height by 1 in. width. The elastic modulus       %
%   is 30x10^6 psi. The frame is subjected to a concenterated load of        %
%   60 lb at the end of the smaller beam and one end of the long member      %
%   is fixed. Use 7 elements to find the deflection of the frame.            %   
%   (see Fig. 8.9.2 for the element discretization)                          %
%                                                                            %
% Variable descriptions                                                      %
%   x and y = global x and y coordiates of each node                         %
%   k = element stiffness matrix                                             %
%   kk = system stiffness matrix                                             %
%   ff = system force vector                                                 %
%   index = a vector containing system dofs associated with each element     %
%   bcdof = a vector containing dofs associated with boundary conditions     %
%   bcval = a vector containing boundary condition values associated with    %
%           the dofs in 'bcdof'                                              %
%----------------------------------------------------------------------------%            

clear
nel=6;           % number of elements
nnel=2;          % number of nodes per element
ndof=3;          % number of dofs per node
nnode=(nnel-1)*nel+1;   % total number of nodes in system
sdof=nnode*ndof; % total system dofs  

x(1)=0; y(1)=0;   % x, y coord. values of node 1 in terms of the global axis
x(2)=0; y(2)=15;  % x, y coord. values of node 2 in terms of the global axis  
x(3)=0; y(3)=30;  % x, y coord. values of node 3 in terms of the global axis
x(4)=0; y(4)=45;  % x, y coord. values of node 4 in terms of the global axis  
x(5)=0; y(5)=60;  % x, y coord. values of node 5 in terms of the global axis
x(6)=10; y(6)=60; % x, y coord. values of node 6 in terms of the global axis  
x(7)=20; y(7)=60; % x, y coord. values of node 7 in terms of the global axis

el=30*10^6;      % elastic modulus
area=2;          % cross-sectional area
xi=2/3;          % moment of inertia of cross-section
rho=1;           % mass density per volume (dummy value for static analysis)
bcdof(1)=1;      % transverse deflection at node 1 is constrained
bcval(1)=0;      % whose described value is 0 
bcdof(2)=2;      % axial displacement at node 1 is constrained
bcval(2)=0;      % whose described value is 0
bcdof(3)=3;      % slope at node 1 is constrained
bcval(3)=0;      % whose described value is 0

ff=zeros(sdof,1);         % initialization of system force vector
kk=zeros(sdof,sdof);      % initialization of system matrix
index=zeros(nel*ndof,1);  % initialization of index vector

ff(20)=-60;     % load applied at node 7 in the negative y direction

for iel=1:nel    % loop for the total number of elements

index=feeldof1(iel,nnel,ndof);  % extract system dofs associated with element

node1=iel;      % starting node number for element 'iel'
node2=iel+1;    % ending node number for element 'iel'

x1=x(node1); y1=y(node1); % x and y coordinate values of 'node1'
x2=x(node2); y2=y(node2); % x and y coordinate values of 'node2'

leng=sqrt((x2-x1)^2+(y2-y1)^2); % length of element 'iel'

if (x2-x1)==0;  % compute the angle between the local and global axes    
	if y2>y1;
   beta=pi/2;
	else
   beta=-pi/2;
   end
else
   beta=atan((y2-y1)/(x2-x1));
end 

k=feframe2(el,xi,leng,area,rho,beta,1); % compute element stiffness matrix

kk=feasmbl1(kk,k,index); % assemble each element matrix into system matrix

end

[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);  % apply the boundary conditions

fsol=kk\ff;   % solve the matrix equation and print

% print both exact and fem solutions


num=1:1:sdof;
store=[num' fsol]


