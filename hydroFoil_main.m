%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Hydroelastic Stability Analysis of an Inverted T Hydrofoil
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialise

close all

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

% Lv = 3; Nv = 12; Lh = 1; Nh = 8; % For debugging only

% Use constant matrial properties (for now)
% E; modulus of elasticity I: second moment of area
tmp = ones(Nv+Nh+Nh-1,1);
%rhoM = 7500;
rhoM = 4000; E = 80e9*tmp; G = 30e9*tmp;
%E=69e9*tmp;
%E=110e9*tmp;
%G=25.5e9*tmp;
%G=42e9*tmp;

sdamp = 100;

% Define the masses
tmp = ones(Nv+Nh+Nh,1); Mmat = zeros(size(tmp)); Imat = zeros(size(tmp));
% Lumped masses at the node points (for translational inertia)
%Mmat = ((0.02*0.02)*1000/4)*tmp;
%Mmat(Nv,1) = Mmat(Nv,1)./2;
% Moment of inertia at the node points (for rotational inertia)
%Imat = 0.001*tmp;

% Number of iterations and eigenmodes to return
% Number of the eigenmode to plot
numits = 300; numeigs = 20; NMplot = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DO NOT MODIFY BELOW HERE
 
% Define the 3D stiffness matrix for an inverted T-shaped hydrofoil

stiffness_main

cntr = 1;
for Uinf = Uinfm
    disp(['Uinf = ' num2str(Uinf)]);
   
    massmat = zeros(size(stiffness));  % Initialize matrices
    Hmass = massmat; Hdamp = Hmass; Hstif = Hmass;

    %------------------------
    % Generate the hydrodynamic integrals for each node
    if Nh < 1
        CLnodes = CLv';
    else
        CLnodes = [CLv';CLh(2:Nh+1)';CLh(2:Nh+1)'];
    end
    
    % Hydrodynamic loading of the vertical section
    Lsmat = zeros(Nv+Nh+Nh,1);
    for i = 1:Nv
        if (i == 1)||(i == Nv)
            Ls = Lv./(2*(Nv-1));
            Lsmat(i,1) = Ls;
        else
            Ls = Lv./(Nv-1);
            Lsmat(i,1) = Ls;
        end
        CL = CLnodes(i);
        % Calculate the centroid of the foil
        [xc,yc] = calcGeom(NACA,HALF_NPANELS,CL);
        % Calculate the hydrodynamic integrals
        [F] = hydroIntegrals(xc,yc,CL,NACA,HALF_NPANELS,rho,Uinf,Ls);
        %F = F.*0;   % Debugging only
        % (For simplicity we will ignore all constant forcing terms)

        % Put these elements into the mass, damping and stiffness matrices
        [Hmass,Hdamp,Hstif] = addHydroIntegralsVertical(Hmass,Hdamp,Hstif,F,i);
    end

    % Hydrodynamic loading of the horizontal section
    % First the mid-section (nvth node)
    if Nh > 0
        Ls = (Lh./Nh);
        CL = CL1;
        [xc,yc] = calcGeom(NACA,HALF_NPANELS,CL);
        [F] = hydroIntegrals(xc,yc,CL,NACA,HALF_NPANELS,rho,Uinf,Ls);
        %F = F.*0;   % Debugging only
        [Hmass,Hdamp,Hstif] = addHydroIntegralsHorizontal(Hmass,Hdamp,Hstif,F,i);

        % Now for the outlying nodes
        for i = (Nv+1):(Nv+Nh+Nh)
            if (i == Nv + Nh)||(i == Nv + Nh + Nh)
                Ls = (Lh./(2*Nh));
                Lsmat(i,1) = Ls;
            else
                Ls = Lh./Nh;
                Lsmat(i,1) = Ls;
            end
            CL = CLnodes(i);
            [xc,yc] = calcGeom(NACA,HALF_NPANELS,CL);
            [F] = hydroIntegrals(xc,yc,CL,NACA,HALF_NPANELS,rho,Uinf,Ls);
            %F = F.*0;   % Debugging only
            [Hmass,Hdamp,Hstif] = addHydroIntegralsHorizontal(Hmass,Hdamp,Hstif,F,i);
        end
    end

    %-------------------------
    % Apply lumped mass parameters
    for i = 1:(Nv+Nh+Nh)
        [xcc,ycc,Ixx,Iyy,JJ,AA] = calcGeom(NACA,HALF_NPANELS,CLnodes(i));
        ind = (i-1)*6 + 1;
        Mmat(i) = Lsmat(i)*rhoM*AA;
        Imat(i) = Lsmat(i)*rhoM*JJ;
        massmat(ind,ind) = massmat(ind,ind) + Mmat(i);
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
    orgmass = massmat;   % Debugging only
    massmat = massmat - Hmass;
    damping = damping + Hdamp;
    stiffness = stiffness + Hstif;

    %------------------------
    % Apply the boundary conditions (all displacements and rotation of node 1 = 0)
    orgstiff = stiffness;
    massmat([1:6],:) = [];massmat(:,[1:6]) = [];
    damping([1:6],:) = [];damping(:,[1:6]) = [];
    stiffness([1:6],:) = [];stiffness(:,[1:6]) = [];

    %------------------------
    % Take out the rows and columns with zero inertia (static condensation)
    nc = size(massmat,1)./6;
    indmat = [];
    invindmat = [];
    for i = 1:size(massmat,1)
        if sum(massmat(:,i)) == 0
            indmat = [indmat,i];
        else
            invindmat = [invindmat,i];
        end
    end

    % sA*theta + sB*x = 0
    % SO: theta = inv(sA)*(-1*sB)*x
    kcc = stiffness(invindmat,invindmat);
    kco = stiffness(invindmat,indmat);
    kcot = stiffness(indmat,invindmat);
    ki = stiffness(indmat,indmat);

    stiff = kcc - kco*(inv(ki)*kcot);
    damp = damping(invindmat,invindmat);
    mass = massmat(invindmat,invindmat);

    %---------------------
    % Reduce to a set of first order differential equations
    %
    % eta_ddot = A*eta_dot - D*eta
    invM = inv(mass);
    A = invM*damp;
    D = invM*stiff;

    % Construct the first order set of equations
    H = [zeros(size(A)) eye(size(A));D A];

    %----------------------
    % Solve the eigenvalue problem
%     [Ve,De] = eig(H);
    opts.maxit = numits;
    opts.disp = 0;
    opts.tol = eps;
    if cntr == 1
        [Ve,De] = eigs(H,numeigs,0+0*i,opts);
    else
        [Ve,De] = eigs(H,numeigs,min(abs(imag(evalsmat([2:size(evalsmat,1)],cntr-1)))),opts);
    end
    
    De = diag(De)
    % Sort the results
    DeR = real(De);DeI=imag(De);
    [DeR,I] = sort(DeR,'descend');DeI = DeI(I,:);
    DeS = complex(DeR,DeI);
    
    if cntr == 1
        evalsmat = zeros(size(DeS,1)+1,size(Uinfm,2));
    end
    evalsmat(:,cntr) = [Uinf;DeS];
    cntr = cntr + 1;

end

figure;
plot(Uinfm,real(evalsmat([2:size(evalsmat,1)],:)),'k.');
grid
xlabel('Foil speed (m/s)')
ylabel('Real part of eigenvalue')

figure;
plot(Uinfm,imag(evalsmat([2:size(evalsmat,1)],:)),'k.');
grid
xlabel('Foil speed (m/s)')
ylabel('Imaginary part of eigenvalue')

% Plot a bode diagram and first mode shape if the length of Uinfm is 1
if (length(Uinfm)==1)
    [G] = getfr(Nv-1,1,Nv-1,1,H,0.1,1000,500,numberNodes,invindmat);
    plotModes3d(0,NMplot,De,Ve,nodeCoordinates,invindmat);
end