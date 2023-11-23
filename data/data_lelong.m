%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lelong hydrofoil (only vertical strut)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho = 1;             % Density of the water (kg/m3)
Uinfm = 1e-6;    % Mean flow velocity (m/s)

NACA = '0015';          % NACA foil shape (4-digit for now)
HALF_NPANELS = 60;      % Half the number of panels for the NACA foil

CL1 = 100e-3;              % Chord length of the foil (at node points)
Lv = 191e-3;              % Length of vertical section
Nv = 20;                % Number of nodes in vertical section (including T vertex)

CL2 = 1e-0; % not needed
Lh = 1e-0;             % Length of horizontal section (2 off)
Nh = 0;                 % Number of nodes in each horizontal section (excluding T vertex)

tmp = ones(Nv+Nh+Nh-1,1);

rhoM = 1420; E = 3.5e9*tmp; G = E/(2*(1+(0.35)));

sdamp = 1;

Mmat = zeros(size(tmp)); Imat = zeros(size(tmp));

dims.Lv=Lv; dims.Lh=Lh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
