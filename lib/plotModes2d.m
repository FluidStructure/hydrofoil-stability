function [] = plotModes2d(dim,dot,N,Ve,cmplx)

%
% Plots the mode shape for a given dimension in the matrix of eigenvectors
% (Ve)
%
% INPUTS:
% dim: 1=x-direction, 2=y-direction, 3=z-direction, 4=torsion
% dot: 0=displacement, 1=velocity
% N = mode number to plot (column of Ve)
% Ve = matrix of eigenvalues
%

ndim = size(Ve,1)/2;

nnods = ndim/4;

mshape = zeros(nnods,1);
for i = 1:nnods
    ind = dot*ndim + ((i-1)*4 + dim);
    mshape(i,1) = Ve(ind,N);
end

if cmplx == 1
    plot(imag([0;mshape]),'rx-')
else
    plot(real([0;mshape]),'rx-')
end
    
    
