%----------------------------------------------------------------------------
% Example 6.5.3                                                              
%   Gauss-Legendre quadrature of a function in 3-dimension
%                                                                            
% Problem description                                                        
%   Integrate f(x,y,z)=1+4x^2y^2-3x^2z^4+y^4z^6 over -1<(x,y,z)<1           
%                                                                            
% Variable descriptions                                                      
%   point3 = integration (or sampling) points                                             
%   weight3 = weighting coefficients                                             
%   nglx = number of integration points along x-axis                                                
%   ngly = number of integration points along y-axis
%   nglz = number of integration points along z-axis
%----------------------------------------------------------------------------%            

clear
nglx=2;           % (2*nglx-1)=2
ngly=3;           % (2*ngly-1)=4
nglz=4;           % (2*nglz-1)=6

[point3,weight3]=feglqd3(nglx,ngly,nglz);  % sampling points & weights

%----------------------------------------------
%  summation for numerical integration
%----------------------------------------------

value=0.0;

for intx=1:nglx
x=point3(intx,1);                  % sampling point in x-axis
wtx=weight3(intx,1);               % weight in x-axis
for inty=1:ngly
y=point3(inty,2);                  % sampling point in y-axis
wty=weight3(inty,2) ;              % weight in y-axis
for intz=1:nglz
z=point3(intz,3);                  % sampling point in z-axis
wtz=weight3(intz,3) ;              % weight in z-axis
func=1+4*x^2*y^2-3*x^2*z^4+y^4*z^6;% evaluate function 
value=value+func*wtx*wty*wtz;
end
end
end

value                    % print the solution

%-----------------------------------------------------------------

