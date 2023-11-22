%
% Plot 3-D contours of the solution to Ex 5.9.2
% See Fig. 5.9.2 for the problem mesh
%
xx=0.0:1.25:5.0;			% x-axis grid
yy=0.0:2.5:10.0;			% y-axis grid

for j=1:5
   for i=1:5
      sol(i,j)=fsol((j-1)*5+i);		% nodal values
   end
end
v=[20 40 60 80];			% values for contour label
c=contour3(xx,yy,sol);			% plot contour
xlabel('x-axis');
ylabel('y-axis');
clabel(c,v);				% label contours

