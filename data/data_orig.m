%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Original data for dubugging

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NACA = '0015';          % NACA foil shape (4-digit for now)
CL1 = 1.3;              % Chord length of the foil (at node points)
CL2 = 0.65;
%CL1 = 0.2;CL2=0.2;
HALF_NPANELS = 60;      % Half the number of panels for the NACA foil

rho = 1000;             % Density of the water (kg/m3)
Uinfm = [0:0.1:40];    % Mean flow velocity (m/s)
Uinfm = 5;

Lv = 3.22;              % Length of vertical section
Nv = 12;                % Number of nodes in vertical section (including T vertex)
Lh = 1.435;             % Length of horizontal section (2 off)
Nh = 8;                 % Number of nodes in each horizontal section (excluding T vertex)
%Nh = 0;
tmp = ones(Nv+Nh+Nh-1,1);

% Lv = 3; Nv = 12; Lh = 1; Nh = 8; % For debugging only

% Use constant matrial properties (for now)
% E; modulus of elasticity I: second moment of area
rhoM = 4000; E = 80e9*tmp; G = 30e9*tmp;
%rhoM = 7500;
%E=69e9*tmp; %E=110e9*tmp;
%G=25.5e9*tmp; %G=42e9*tmp;

sdamp = 100;

% Define the masses
Mmat = zeros(size(tmp)); Imat = zeros(size(tmp));
% Lumped masses at the node points (for translational inertia)
%Mmat = ((0.02*0.02)*1000/4)*tmp;
%Mmat(Nv,1) = Mmat(Nv,1)./2;
% Moment of inertia at the node points (for rotational inertia)
%Imat = 0.001*tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
