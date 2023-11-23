function [] = plotModes3d(dot,N,De,Ve,nodeCoordinates,invindmat,dims)

%
% Plots the mode shape for all translational directions in a 3dplot for
% a given dimension in the matrix of eigenvectors (Ve) OVER ONE COMPLETE CYCLE
% OF OSCILLATION
%
% INPUTS:
% dot: 0=displacement, 1=velocity
% N = mode number to plot (column of Ve)
% De = matrix of eigenvalues 
% Ve = matrix of eigenvectors
%
% Scaling factor (maximum amplitude of plot)
SCF = 0.3;

% save the individual nodecoordinates
nodeCoordinates(1,:) = [];
xx = nodeCoordinates(:,1);
yy = nodeCoordinates(:,2);
zz = nodeCoordinates(:,3);

% Get the number of dimensions and the number of nodes
ndim = size(Ve,1)/2;
nnods = size(nodeCoordinates,1);

% Re-fill the eigenvector with zero-inertia DOFs (originally removed)
Vet = Ve([1:ndim],N);
Ven = zeros(nnods*6,1);
Ven(invindmat) = Vet;

% Intialise the figure window
fig1 = figure;

% Plot the base node coordinate positions
%subplot(2,1,1);
plot3(xx,yy,zz,'b.');
%grid;rotate3d;

% Get the eigenvalue
ev = De(N);
% SUPPRESS GROWTH/DECAY COMPONENT!!
ev = complex(0,imag(ev));

% Get the eigenvectors
dx = zeros(nnods,1);
dy = dx;
dz = dx;
for i = 1:nnods
    ind = dot*ndim + ((i-1)*6 + 1);
    dx(i,1) = Ven(ind,1);
    dy(i,1) = Ven(ind+1,1);
    dz(i,1) = Ven(ind+2,1);
end

% Get the time period of oscillation
T = 2*pi./imag(ev);

% Plot the position of the wall over one time period of oscillation
%subplot(2,1,2)
tmat = [0:(T/81):T];
x = zeros(length(dx),length(tmat));
y = x;z = x;
c = 1;
for tt = tmat
    x(:,c) = real(dx.*exp(ev.*tt));
    y(:,c) = real(dy.*exp(ev.*tt));
    z(:,c) = real(dz.*exp(ev.*tt));
    c = c + 1;
end
% Scale the results
mdef = max(max([x y z]));
x = (x./mdef)*SCF;
y = (y./mdef)*SCF;
z = (z./mdef)*SCF;
% Plot the results
c = 1;
hold on;
for tt = tmat
    if c == 1
        plot3(x(:,c)+xx,y(:,c)+yy,z(:,c)+zz,'ko')
    else
        plot3(x(:,c)+xx,y(:,c)+yy,z(:,c)+zz,'k.')
    end
    c = c + 1;
end
hold off;
axis([-dims.Lh dims.Lh -dims.Lh dims.Lh -dims.Lv 0])
%axis([-2 2 -2 2 -4 0])
grid;rotate3d;

