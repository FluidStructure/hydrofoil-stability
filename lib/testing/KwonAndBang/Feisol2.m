function [shape,dhdr]=kwisol2(rvalue)
%--------------------------------------------------------
%Purpose
%	Compute isoparametric 2-node shape functions
%	and their derivatives at the selected
%       point in terms of the natural coordinate.
%   
%Synopsis:
%	[shape,dhdr]=kwisol2(rvalue)
%
%Variable Description:
%	shape - shape functions for the linear element
%       dhdr - derivatives of shape functions
%       rvalue - r coordinate valute of the selected point   
%   
%Notes:
%	1st node at rvalue=-1
%       2nd node at rvalue=1
%--------------------------------------------------------

shape(1)=0.5*(1.0-rvalue);  % first shape function
shape(2)=0.5*(1.0+rvalue);  % second shape function

dhdr(1)=-0.5;  % derivative of the first shape function
dhdr(2)=0.5; % derivative of the second shape function


