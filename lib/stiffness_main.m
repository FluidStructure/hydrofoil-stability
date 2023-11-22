%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the 3D stiffness matrix for an inverted T-shaped hydrofoil

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generation of coordinates and connectivities

vsec = linspace(0,-1*Lv,Nv)';
hsec = linspace(0,Lh,(Nh+1))';hsec(1) = [];
% Coordinates: Z = vertical, X = horizontal (side<->side), Y = for<->aft
nodeCoordinates=[zeros(size(vsec)) zeros(size(vsec)) vsec];
%nodeCoordinates=[-1*vsec zeros(size(vsec)) zeros(size(vsec))];      % For debugging only
nodeCoordinates=[nodeCoordinates;[hsec zeros(size(hsec)) min(vsec).*ones(size(hsec))]];
nodeCoordinates=[nodeCoordinates;[-1*hsec zeros(size(hsec)) min(vsec).*ones(size(hsec))]];

xx=nodeCoordinates(:,1); yy=nodeCoordinates(:,2);
zz=nodeCoordinates(:,3); numberNodes=size(nodeCoordinates,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the foil chord length

CLv = linspace(CL1,CL2,Nv); CLh = linspace(CL1,CL2,Nh+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define connectivity

NEv = Nv - 1;   % Number of elements in vertical section
NEh = Nh;       % Number of elements in each horizontal section
ec = 1;

% Vertical section
elementNodes = zeros((NEv + 2*NEh),2);
CLe = zeros((NEv + 2*NEh),1);
for i = 1:NEv
    elementNodes(ec,[1 2])=[ec ec+1];
    CLe(ec,1) = (CLv(ec) + CLv(ec+1))/2;
    ec = ec + 1;
end
ecv = ec - 1;
if Nh > 0
    % Horizontal element #1
    elementNodes(ec,[1 2])=[ec+1 NEv+1];
    CLe(ec) = (CLh(1) + CLh(2))/2;
    ec = ec + 1;
    for i = 1:(NEh-1)
        elementNodes(ec,[1 2])=[ec+1 ec];
        CLe(ec) = (CLh(ec - ecv) + CLh(ec - ecv + 1))/2;
        ec = ec + 1;
    end
    ecv = ec - 1;
    % Horizontal element #2
    elementNodes(ec,[1 2])=[ec+1 NEv+1];
    CLe(ec) = (CLh(1) + CLh(2))/2;
    ec = ec + 1;
    for i = 1:(NEh-1)
        elementNodes(ec,[1 2])=[ec+1 ec];
        CLe(ec) = (CLh(ec - ecv) + CLh(ec - ecv + 1))/2;
        ec = ec + 1;
    end
end

numberElements=size(elementNodes,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the centroid & moments of inertia (xx, yy and polar) for the elements

Amat=zeros(size(E));       % A=(0.1*0.1)*tmp;
Iy=zeros(size(E));      % Iy=(0.1*(0.1^3)/12)*tmp;
Iz=zeros(size(E));      % Iz=(0.1*(0.1^3)/12)*tmp;
J=zeros(size(E));       % J=5e-5*tmp;
xce=zeros(size(E));
yce=zeros(size(E));
for i = 1:numberElements
    [xcc,ycc,Ixx,Iyy,JJ,AA] = calcGeom(NACA,HALF_NPANELS,CLe(i));
    xce(i) = xcc; yce(i) = ycc;
    Iy(i) = Ixx;  Iz(i) = Iyy;
    J(i) = JJ;  Amat(i) = AA;
end

% Amat = Amat + (0.02*0.02);
% Iy = Iy + 100*0.02*(0.02^3)/12;
% Iz = Iz + 0.02*(0.02^3)/12;
% J = J + 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate the stiffness, mass and damping matix

GDof=6*numberNodes;
% E = ones(size(E));Amat=E;Iz=ones(size(E));Iy=100*ones(size(E));G=E;J=100*ones(size(E));
[stiffness]=formStiffness3Dframe(GDof,numberElements,elementNodes,...
                                 nodeCoordinates,E,Amat,Iz,Iy,G,J);
stiffness = stiffness*-1;
damping   = -1*sdamp*eye(size(stiffness));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
