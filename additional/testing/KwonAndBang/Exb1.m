%
% generate uniform 2-D rectangular mesh for a rectangular
% shape of domain
%
clear
xnode=5;  	% number of nodes in x-axis
ynode=9;	% number of nodes in y-axis

xzero=0.0;	% x-coordniate of the bottom left corner
yzero=0.0;	% y-coordinate of the bottom left corner

xlength=2.0;	% size of domain in x-axis
ylength=4.0; % size of domain in y-axis

xnel=xnode-1;	% number of elements in x-axis
ynel=ynode-1;	% number of elements in y-axis

nel=xnel*ynel;	% total number of elements
nnode=xnode*ynode;	% total number of nodes

delx=xlength/xnel;	% element size in x-axis
dely=ylength/ynel;	% element size in y-axis

nodes=zeros(nel,4);	
gcoord=zeros(nnode,2);

% generate nodal coordinate values

for iy=1:ynode
for ix=1:xnode
gcoord((iy-1)*xnode+ix,1)=xzero+(ix-1)*delx;
gcoord((iy-1)*xnode+ix,2)=yzero+(iy-1)*dely;
end
end

% generate nodal connectivity in ccw direction
 
for iiy=1:ynel
for iix=1:xnel
nodes((iiy-1)*xnel+iix,1)=(iiy-1)*xnode+iix;
nodes((iiy-1)*xnel+iix,2)=(iiy-1)*xnode+iix+1;
nodes((iiy-1)*xnel+iix,3)=iiy*xnode+iix+1;
nodes((iiy-1)*xnel+iix,4)=iiy*xnode+iix;
end
end


