%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lelong hydrofoil (only vertical strut)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho = 1000;             % Density of the water (kg/m3)
Uinfm = 0;    % Mean flow velocity (m/s)

NACA = '0006';
if rho==1;NACA = '0002';end
HALF_NPANELS = 60;      % Half the number of panels for the NACA foil

CL1 = 100e-3;              % Chord length of the foil (at node points)
CL2 = 100e-3;

Lv = 191e-3;              % Length of vertical section
Nv = 20;                % Number of nodes in vertical section (including T vertex)

Lh = 0;             % Length of horizontal section (2 off)
Nh = 0;                 % Number of nodes in each horizontal section (excluding T vertex)

tmp = ones(Nv+Nh+Nh-1,1);

rhoM = 1420; E = 3.5e9*tmp; G = E/(2*(1+(0.35)));

sdamp = 0;

Mmat = zeros(size(tmp)); Imat = zeros(size(tmp));

dims.Lv=Lv; dims.Lh=Lh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
