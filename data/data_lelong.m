%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lelong hydrofoil (only vertical strut)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NACA = '0015';          % NACA foil shape (4-digit for now)
CL1 = 1.3;              % Chord length of the foil (at node points)
CL2 = 0.65;
HALF_NPANELS = 60;      % Half the number of panels for the NACA foil

rho = 1000;             % Density of the water (kg/m3)
Uinfm = [0:0.1:40];    % Mean flow velocity (m/s)
Uinfm = 5;

Lv = 3.22;              % Length of vertical section
Nv = 12;                % Number of nodes in vertical section (including T vertex)
Lh = 1.435;             % Length of horizontal section (2 off)
Nh = 0;                 % Number of nodes in each horizontal section (excluding T vertex)

tmp = ones(Nv+Nh+Nh-1,1);

rhoM = 4000; E = 80e9*tmp; G = 30e9*tmp;

sdamp = 0;

Mmat = zeros(size(tmp)); Imat = zeros(size(tmp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
