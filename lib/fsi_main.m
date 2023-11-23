%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% couple fluid and structure & solve system

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply lumped mass parameters

massmat = zeros(size(stiffness));  % Initialize matrices
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add hydrodynamic integrals to the equations

%Hdamp = Hdamp*0;    % Debugging only
%orgmass   = massmat;   % Debugging only
massmat   = massmat   - Hmass;
damping   = damping   + Hdamp;
stiffness = stiffness + Hstif;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply the boundary conditions (all displacements and rotation of node 1 = 0)

%orgstiff = stiffness;   % Debugging only
massmat  ([1:6],:) = []; massmat  (:,[1:6]) = [];
damping  ([1:6],:) = []; damping  (:,[1:6]) = [];
stiffness([1:6],:) = []; stiffness(:,[1:6]) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Take out the rows and columns with zero inertia (static condensation)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reduce to a set of first order differential equations

% eta_ddot = A*eta_dot - D*eta
invM = inv(mass); A = invM*damp;  D = invM*stiff;
%A = damp/mass;  D = stiff/mass; % same soln, but mode plots different!

% Construct the first order set of equations
H = [zeros(size(A)) eye(size(A));D A];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve the eigenvalue problem

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

