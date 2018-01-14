%----------------------------------------------------------------------------
% Example 6.5.1                                                              
%   Gauss-Legendre quadrature of a function in 1-dimension
%                                                                            
% Problem description                                                        
%   Integrate f(x)=1+x^2-3x^3+4x^5 between x=-1 and x=1           
%                                                                            
% Variable descriptions                                                      
%   point1 = integration (or sampling) points                                             
%   weight1 = weighting coefficients                                             
%   ngl = number of integration points                                                 
%----------------------------------------------------------------------------%            

clear
ngl=3;           % (2*ngl-1)=5

[point1,weight1]=feglqd1(ngl);  % extract integration points ans weights

%----------------------------------------------
%  summation for numerical integration
%----------------------------------------------

value=0.0;

for int=1:ngl
x=point1(int);
wt=weight1(int);
func=1+x^2-3*x^3+4*x^5;    % evaluate function at the integration point
value=value+func*wt;
end

value                    % print the solution

%-----------------------------------------------------------------

