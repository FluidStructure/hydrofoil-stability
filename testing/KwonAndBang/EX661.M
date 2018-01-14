%-----------------------------------------------------------------------------
% Example 6.6.1                                                              
%   Compute element matrix for two-dimensional Laplace equation
%                                                                            
% Problem description                                                        
%   Determine the element matrix for Laplace equation using 
%   isoparametric four-node quadrilateral element and Gauss-Legendre          
%   quadrature for a single element shown in Fig. 6.2.4.
%                                                                            
% Variable descriptions                                                      
%   k - element matrix
%   point2 - integration (or sampling) points                                             
%   weight2 - weighting coefficients                                             
%   nglx - number of integration points along x-axis                                                
%   ngly - number of integration points along y-axis
%   xcoord - x coordinate values of nodes
%   ycoord - y coordinate values of nodes
%   jacob2 - jacobian matrix
%   shape - four-node quadrilateral shape functions
%   dhdr - derivatives of shape functions w.r.t. natural coord. r
%   dhds - derivatives of shape functions w.r.t. natural coord. s
%   dhdx - derivatives of shape functions w.r.t. physical coord. x
%   dhdy - derivatives of shape functions w.r.t. physical coord. y
%-----------------------------------------------------------------------------            

clear
nnel=4;                               % number of nodes per element
ndof=1;                               % degrees of freedom per node
edof=nnel*ndof;                       % degrees of freedom per element

nglx=2; ngly=2;                       % use 2x2 integration rule

xcoord=[-1 1 1 -1];                   % x coordinate values
ycoord=[-0.75 -0.75 1.25 0.25];       % y coordinate values

[point2,weight2]=feglqd2(nglx,ngly);  % sampling points & weights

%--------------------------------
%  numerical integration
%--------------------------------

k=zeros(edof,edof);                   % initialization to zero

for intx=1:nglx
x=point2(intx,1);                  % sampling point in x-axis
wtx=weight2(intx,1);               % weight in x-axis
for inty=1:ngly
y=point2(inty,2);                  % sampling point in y-axis
wty=weight2(inty,2) ;              % weight in y-axis

[shape,dhdr,dhds]=feisoq4(x,y); % compute shape functions and
                                      % derivatives at sampling point

jacob2=fejacob2(nnel,dhdr,dhds,xcoord,ycoord);  % compute Jacobian

detjacob=det(jacob2);                 % determinant of Jacobian
invjacob=inv(jacob2);                 % inverse of Jacobian matrix

[dhdx,dhdy]=federiv2(nnel,dhdr,dhds,invjacob); % derivatives w.r.t.
                                                   % physical coordinate

%------------------------------
%  element matrix loop
%------------------------------

for i=1:edof
for j=1:edof
k(i,j)=k(i,j)+(dhdx(i)*dhdx(j)+dhdy(i)*dhdy(j))*wtx*wty*detjacob;
end
end

end
end

k                        % print the element matrix
%-----------------------------------------------------------------

