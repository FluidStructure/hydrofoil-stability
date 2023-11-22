function [num, den] = festotf(A,B,C,D,iu)
%---------------------------------------------------------------------------------------------
%  Purpose:
%     The funcstion subroutine festotf.m converts a state space form of system
%     into a transfer function form.
%
%     For a given system
%      
%               xdot = Ax+Bu,
%               y =Cx+Du
% 
%     The transfer function becomes
%
%			N(s)          -1
%		H(s) = -------- = C(sI-A) B + D
%			D(s)
%  Synopsis:
%     [num,den]=festotf(A,B,C,D,iu)       
%
%  Variable Description:       
%     Input parameters  - System matries [A,B,C,D]
%                         iu - Index for control input(iu-th input)
%	
%     Output  parameters - D(s)  -  Vector of coefficients of the denominator polynomial
%                          N(s) - Vector of coefficients of the numerator polynomials
%	 	                - There are same number of rows in N(s) as the number of output
%-----------------------------------------------------------------------------------------------


den = poly(A);						%   Determine denominator polynomial

B = B(:,iu);			                        %  Select the corresponding column 
				
D = D(:,iu);
[m,n] = size(C);
num = ones(m, n+1);
for i=1:m
	num(i,:) = poly(A-B*C(i,:)) + (D(i) - 1) * den;
end
