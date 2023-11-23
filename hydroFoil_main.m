%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Hydroelastic Stability Analysis of an Inverted T Hydrofoil
% M. Pitman 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialise

clear; close all

if ~exist('getfr.m','file')
    addpath(genpath('lib')); % addpath(genpath('data'));
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data

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

% Number of iterations and eigenmodes to return & the eigenmode to plot
numits = 300; numeigs = 20; NMplot = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% main code

cntr = 1;
for Uinf = Uinfm
    disp(['Uinf = ' num2str(Uinf)]);

    stiffness_main % 3D stiffness matrix for an inverted T-shaped hydrofoil

    massmat = zeros(size(stiffness));  % Initialize matrices
    Hmass = massmat; Hdamp = Hmass; Hstif = Hmass;

    hydro_main % Generate the hydrodynamic integrals for each node

    %-------------------------
    % Apply lumped mass parameters
    for i = 1:(Nv+Nh+Nh)
        [xcc,ycc,Ixx,Iyy,JJ,AA] = calcGeom(NACA,HALF_NPANELS,CLnodes(i));
        ind = (i-1)*6 + 1;
        Mmat(i) = Lsmat(i)*rhoM*AA;   Imat(i) = Lsmat(i)*rhoM*JJ;
        massmat(ind  ,ind  ) = massmat(ind  ,ind  ) + Mmat(i);
        massmat(ind+1,ind+1) = massmat(ind+1,ind+1) + Mmat(i);
        massmat(ind+2,ind+2) = massmat(ind+2,ind+2) + Mmat(i);
        if i <= Nv
            massmat(ind+5,ind+5) = massmat(ind+5,ind+5) + Imat(i);
        else
            massmat(ind+3,ind+3) = massmat(ind+3,ind+3) + Imat(i);
        end
    end

    %------------------------
    % Add hydrodynamic integrals to the equations
    %Hdamp = Hdamp*0;    % Debugging only
    %orgmass   = massmat;   % Debugging only
    massmat   = massmat   - Hmass;
    damping   = damping   + Hdamp;
    stiffness = stiffness + Hstif;

    %------------------------
    % Apply the boundary conditions (all displacements and rotation of node 1 = 0)
    %orgstiff = stiffness;   % Debugging only
    massmat  ([1:6],:) = []; massmat  (:,[1:6]) = [];
    damping  ([1:6],:) = []; damping  (:,[1:6]) = [];
    stiffness([1:6],:) = []; stiffness(:,[1:6]) = [];

    %------------------------
    % Take out the rows and columns with zero inertia (static condensation)
    nc = size(massmat,1)./6; indmat = []; invindmat = [];
    for i = 1:size(massmat,1)
        if sum(massmat(:,i)) == 0
            indmat = [indmat,i];
        else
            invindmat = [invindmat,i];
        end
    end

    % sA*theta + sB*x = 0... SO: theta = inv(sA)*(-1*sB)*x
    kcc  = stiffness(invindmat,invindmat);
    kco  = stiffness(invindmat,indmat   );
    kcot = stiffness(indmat   ,invindmat);
    ki   = stiffness(indmat   ,indmat   );

    stiff = kcc - kco*(inv(ki)*kcot);
    damp  = damping(invindmat,invindmat);
    mass  = massmat(invindmat,invindmat);

    %---------------------
    % Reduce to a set of first order differential equations
    %
    % eta_ddot = A*eta_dot - D*eta
    invM = inv(mass); A = invM*damp;  D = invM*stiff;
    %A = damp/mass;  D = stiff/mass; % same soln, but mode plots different!

    % Construct the first order set of equations
    H = [zeros(size(A)) eye(size(A));D A];

    %----------------------
    % Solve the eigenvalue problem

    opts.maxit = numits; opts.disp = 0; opts.tol = eps;
    if cntr == 1
        [Ve,De] = eigs(H,numeigs,0+0*i,opts);
    else
        [Ve,De] = eigs(H,numeigs,min(abs(imag(evalsmat([2:size(evalsmat,1)],cntr-1)))),opts);
    end
    De = diag(De)

    % Sort the results
    DeR = real(De); DeI = imag(De);
    [DeR,I] = sort(DeR,'descend'); DeI = DeI(I,:);
    DeS = complex(DeR,DeI);
    
    if cntr == 1
        evalsmat = zeros(size(DeS,1)+1,size(Uinfm,2));
    end
    evalsmat(:,cntr) = [Uinf;DeS];
    cntr = cntr + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot

figure;
plot(Uinfm,real(evalsmat([2:size(evalsmat,1)],:)),'k.'); grid
xlabel('Foil speed (m/s)'); ylabel('Real part of eigenvalue')

figure;
plot(Uinfm,imag(evalsmat([2:size(evalsmat,1)],:)),'k.'); grid
xlabel('Foil speed (m/s)'); ylabel('Imaginary part of eigenvalue')

% Plot a bode diagram and first mode shape if the length of Uinfm is 1
if (length(Uinfm)==1)
    [G] = getfr(Nv-1,1,Nv-1,1,H,0.1,1000,500,numberNodes,invindmat);
    plotModes3d(0,NMplot,De,Ve,nodeCoordinates,invindmat);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
