%----------------------------------------------------------------------------
% Example 6.5.2                                                              
%   Gauss-Legendre quadrature of a function in 2-dimension
%                                                                            
% Problem description                                                        
%   Integrate f(x,y)=1+4xy-3x^2y^2+x^4y^6 over -1<x<1 and -1<y<1          
%                                                                            
% Variable descriptions                                                      
%   point2 = integration (or sampling) points                                             
%   weight2 = weighting coefficients                                             
%   nglx = number of integration points along x-axis                                                
%   ngly = number of integration points along y-axis
%----------------------------------------------------------------------------%            

clear
nglx=3;           % (2*nglx-1)=4
ngly=4;           % (2*ngly-1)=6

[point2,weight2]=feglqd2(nglx,ngly);  % integration points and weights

%----------------------------------------------
%  summation for numerical integration
%----------------------------------------------

value=0.0;

for intx=1:nglx
x=point2(intx,1);                  % sampling point in x-axis
wtx=weight2(intx,1);               % weight in x-axis
for inty=1:ngly
y=point2(inty,2);                  % sampling point in y-axis
wty=weight2(inty,2) ;              % weight in y-axis
func=1+4*x*y-3*x^2*y^2+x^4*y^6;    % evaluate function at sampling points
value=value+func*wtx*wty;
end
end

value                    % print the solution

%-----------------------------------------------------------------

