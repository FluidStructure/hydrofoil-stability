%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Curtin aluminium plate (only vertical strut)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho = 1e0;             % Density of the water (kg/m3)
Uinfm = [0:0.5:20];    % Mean flow velocity (m/s)
%Uinfm = 0;    % Mean flow velocity (m/s)

NACA = '6407';          % NACA foil shape (4-digit for now)
HALF_NPANELS = 60;      % Half the number of panels for the NACA foil

CL1 = 50e-3;              % Chord length of the foil (at node points)
CL2 = 50e-3;

Lv = 440e-3;              % Length of vertical section
Nv = 20;                % Number of nodes in vertical section (including T vertex)

Lh = 0;             % Length of horizontal section (2 off)
Nh = 0;                 % Number of nodes in each horizontal section (excluding T vertex)

tmp = ones(Nv+Nh+Nh-1,1);

rhoM = 2720; E = 68.9e9*tmp; G = 3.75e9*tmp;

sdamp = 0;

Mmat = zeros(size(tmp)); Imat = zeros(size(tmp));

dims.Lv=Lv; dims.Lh=Lh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
