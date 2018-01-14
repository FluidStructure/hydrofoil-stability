%----------------------------------------------------------------------------
% Example 12.9.1                                                              
% to solve the ordinary differential equation given as            
%   u'' = 1,  0 < x < 1                                                                        
%   u(0) = 0  and  u(1) = 0
% and to compute the mean-square norm for different element sizes.
%
% Variable descriptions                                                      
%   k = element matrix                                             
%   f = element vector
%   kk = system matrix                                             
%   ff = system vector                                                 
%   index = a vector containing system dofs associated with each element     
%   bcdof = a vector containing dofs associated with boundary conditions     
%   bcval = a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'                                              
%----------------------------------------------------------------------------            

%------------------------------------
%  input data for control parameters
%------------------------------------

clear
for int=1:7      % loop for the different element numbers

nel=2^int;                  % number of elements
nnel=2;                 % number of nodes per element
ndof=1;                 % number of dofs per node
nnode=nel+1;                % total number of nodes in system
sdof=nnode*ndof;        % total system dofs  

%-----------------------------------------
%  input data for nodal coordinate values
%-----------------------------------------

elemsize(int)=1.0/nel;
for i=1:nnode
gcoord(i)=elemsize(int)*(i-1);
end

%-----------------------------------------------------
%  input data for nodal connectivity for each element
%-----------------------------------------------------

for i=1:nel
nodes(i,1)=i;  
nodes(i,2)=i+1; 
end

%-----------------------------------------
%  input data for coefficients of the ODE
%-----------------------------------------

acoef=1;                % coefficient 'a' of the diff eqn
bcoef=0;               % coefficient 'b' of the diff eqn
ccoef=0;                % coefficient 'c' of the diff eqn    

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------

bcdof(1)=1;             % first node is constrained
bcval(1)=0;             % whose described value is 0 
bcdof(2)=nnode;             % 4th node is constrained
bcval(2)=0;             % whose described value is 0

%-----------------------------------------
%  initialization of matrices and vectors
%-----------------------------------------

ff=zeros(sdof,1);       % initialization of system force vector
kk=zeros(sdof,sdof);    % initialization of system matrix
index=zeros(nnel*ndof,1);  % initialization of index vector

%-----------------------------------------------------------------
%  computation of element matrices and vectors and their assembly
%-----------------------------------------------------------------

for iel=1:nel           % loop for the total number of elements

nl=nodes(iel,1); nr=nodes(iel,2); % extract nodes for (iel)-th element
xl=gcoord(nl); xr=gcoord(nr);% extract nodal coord values for the element
eleng=xr-xl;            % element length
index=feeldof1(iel,nnel,ndof);% extract system dofs associated with element

k=feode2l(acoef,bcoef,ccoef,eleng); % compute element matrix
f=fef1l(xl,xr);                     % compute element vector
[kk,ff]=feasmbl2(kk,ff,k,f,index);  % assemble element matrices and vectors

end		% end of loop for elments

%-----------------------------
%  apply boundary conditions
%-----------------------------

[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);

%----------------------------
%  solve the matrix equation
%----------------------------

fsol=kk\ff;   

%-----------------------------
%  compute mean-square error
%-----------------------------

error(int)=0.0;

for iel=1:nel

nl=nodes(iel,1); nr=nodes(iel,2); % extract nodes for (iel)-th element
xl=gcoord(nl); xr=gcoord(nr);% extract nodal coord values for the element
eleng=xr-xl;            % element length
soll=fsol(nl); solr=fsol(nr); % extract fem solution

co1=(soll-solr)/eleng-0.5;
co2=(xl*solr-xr*soll)/eleng;

error(int)=error(int)+(xr^5-xl^5)/20.0+co1*(xr^4-xl^4)/4.0 ...
 +(co1^2+co2)*(xr^3-xl^3)/3.0+co1*co2*(xr^2-xl^2)+co2^2*eleng;

end

error(int)=sqrt(error(int));

end % end of loop of different element numbers

%------------------------------------------------------------
%  log-log plot of the mean-square error vs element size
%------------------------------------------------------------

loglog(elemsize,error);
xlabel('Element Size')
ylabel('Mean-Square Error')
