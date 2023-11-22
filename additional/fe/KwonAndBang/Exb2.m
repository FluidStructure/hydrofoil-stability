
% x and y coordinates

xcoord=gcoord(:,1);
ycoord=gcoord(:,2);

% extract coordinates for each element

for i=1:nel	
for j=1:4
x(j)=xcoord(nodes(i,j));
y(j)=ycoord(nodes(i,j));
end; %j loop
	
xvec=[x(1),x(2),x(3),x(4),x(1)];	
yvec=[y(1),y(2),y(3),y(4),y(1)];	
plot(xvec,yvec);	
hold on;
	
% place element number
		
midx=mean(xvec(1:4));
midy=mean(yvec(1:4));		
text(midx,midy,num2str(i));

end; % i loop

xlabel('x axis');
ylabel('y axis');
title(num2str(i));

% put node numbers

for jj=1:nnode
text(gcoord(jj,1),gcoord(jj,2),num2str(jj));
end;

