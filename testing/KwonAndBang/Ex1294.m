%----------------------------------------------------------------------------
% Example 12.9.4                                                              
% to solve the ordinary differential equation given as            
%   u'' - alpa*u*u'= 0,  0 < x < 1                                                                        
%   u(0) = 0  and  u(1) = 1
% using classical linearization (u')*u.
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

nel=100;                  % number of elements
nnel=2;                 % number of nodes per element
ndof=1;                 % number of dofs per node
nnode=nel+1;                % total number of nodes in system
sdof=nnode*ndof;        % total system dofs  
alpa=100;      % coefficient of the nonlinear term
toler=0.0001;   % error tolerance to terminate iterations

%-----------------------------------------
%  input data for nodal coordinate values
%-----------------------------------------

elemsize=1.0/nel;
for i=1:nnode
gcoord(i)=elemsize*(i-1);
end

%-----------------------------------------------------
%  input data for nodal connectivity for each element
%-----------------------------------------------------

for i=1:nel
nodes(i,1)=i;  
nodes(i,2)=i+1; 
end

%-------------------------------------
%  input data for boundary conditions
%-------------------------------------

bcdof(1)=1;             % first node is constrained
bcval(1)=0;             % whose described value is 0 
bcdof(2)=nnode;             % 4th node is constrained
bcval(2)=1;             % whose described value is 0

%-----------------------------------------
%  loop for iteration
%-----------------------------------------

error=1;   % error is set to 1 arbitrarily
solold=gcoord'; % assume a linear function initially

it=0;
while error > toler
it=it+1;  % iteration counter

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
soll=solold(nl); solr=solold(nr); % extract old solutions

%---------------------------------------------------
%  element matrix and vector for quasilinearization
%---------------------------------------------------

k(1,1)=1/eleng+alpa*(solr-soll)/3; % element matrix
k(1,2)=-1/eleng+alpa*(solr-soll)/6;
k(2,1)=-1/eleng+alpa*(solr-soll)/6;
k(2,2)=1/eleng+alpa*(solr-soll)/3;

f(1)=0;  % element vector
f(2)=0;

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

%----------------------------
%  check the error
%----------------------------

error=0;
for i=2:(nnode-1)
error=error+(fsol(i)-solold(i))^2/fsol(i)^2;
end

error=sqrt(error);

solold=fsol;    %  assigne previous solution

end  % end of iteration loop

%--------------------------------
%  plot of the solution
%--------------------------------

No_of_Iteration=it % print number of iteration for convergency

plot(gcoord',fsol);
xlabel('x-axis')
ylabel('Solution')
 
%-------------------------------

