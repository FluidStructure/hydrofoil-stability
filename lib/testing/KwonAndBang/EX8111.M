%--------------------------------------------------------------------------
% Example 8.11.1                                                           
% to find the transient response of a cantilever beam with a tip load
%                                                                          
% Problem description                                                        
%   Find the transient response of a cantilever beam whose length is 1 m.      
%   long. The beam has the cross-section of 0.02 m by 0.02 m and the mass   
%   density is 1000 Kg/m^3. The elastic and shear modulus is 100 GPAa and         
%   40 GPa, respectively. Use 4 elements.                                
%                                                                            
% Variable descriptions                                                      
%   k = element stiffness matrix                                             
%   kk = system stiffness matrix                                             
%   m = element mass matrix                                                  
%   mm = system mass matrix                                                  
%   index = a vector containing system dofs associated with each element     
%   bcdof = a vector containing dofs associated with boundary conditions     
%------------------------------------------------------------------------ 

clear
nel=4;              % number of elements
nnel=2;             % number of nodes per element
ndof=2;             % number of dofs per node
nnode=(nnel-1)*nel+1;   % total number of nodes in system
sdof=nnode*ndof;    % total system dofs  

el=100*10^9;        % elastic modulus
tleng=1;            % total beam length
leng=tleng/nel;     % same size of beam elements
xi=0.02^4/12;       % height (or thickness) of the beam
area=0.004;	    % cross-sectional area of the beam
rho=1000;           % mass density of the beam
ipt=1;		    % option flag for mass matrix

dt=0.0001;	    % time step size
ti=0;               % initial time
tf=0.2;             % final time
nt=fix((tf-ti)/dt); % number of time steps

nbc=2;              % number of constraints
bcdof(1)=1;         % transverse displ. at node 1 is constrained
bcdof(2)=2;         % slope at node 1 is constrained

kk=zeros(sdof,sdof);  % initialization of system stiffness matrix
mm=zeros(sdof,sdof);  % initialization of system mass matrix
force=zeros(sdof,1);  % initialization of force vector
index=zeros(nel*ndof,1);  % initialization of index vector
acc=zeros(sdof,nt);   % initialization of acceleartion matrix
vel=zeros(sdof,nt);   % initialization of velocity matrix
disp=zeros(sdof,nt); % initialization of displ. matrix

vel(:,1)=zeros(sdof,1);         % initial zero velocity
disp(:,1)=zeros(sdof,1);        % initial zero displacement
force(9)=100;     % tip load of 100 

for iel=1:nel    % loop for the total number of elements

index=feeldof1(iel,nnel,ndof); % extract system dofs associated with element

[k,m]=febeam1(el,xi,leng,area,rho,ipt); % compute element stiffness matrix

kk=feasmbl1(kk,k,index); % assemble each element matrix into system matrix

mm=feasmbl1(mm,m,index); % assemble each element matrix into system matrix

end

mminv=inv(mm);        % invert the mass matrix

% central difference scheme for time integration

for it=1:nt

acc(:,it)=mminv*(force-kk*disp(:,it));

  for i=1:nbc
  ibc=bcdof(i);  
  acc(ibc,it)=0;
  end

vel(:,it+1)=vel(:,it)+acc(:,it)*dt;
disp(:,it+1)=disp(:,it)+vel(:,it+1)*dt;

end

acc(:,nt+1)=mminv*(force-kk*disp(:,nt+1));

time=0:dt:nt*dt;
plot(time,disp(9,:))
xlabel('Time(seconds)')
ylabel('Tip displ. (m)')










































