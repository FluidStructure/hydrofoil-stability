%
% Constructs the linear equations for inverted "T" hydrofoil using lumped
% mass and linear hydrodynamic stability analysis
%

% Four mass points like this:
%||---x---x---x---x

nnods = 5;

% Input values of EI between lumped mass points
%EImat = [1;1;1];
EImat = ones(nnods,1)*1;
% Input length between
%Lmat = [0.5;0.5;0.5];
Lmat = ones(nnods,1)./nnods;
% Input masses
%Mmat = [1;1;0.5];
Mmat = ones(nnods,1);Mmat(nnods,1) = Mmat(nnods,1)./2;Mmat = Mmat*3./nnods;

%----------------------
% Generate the stiffness matrices for the vertical section in the
% transverse direction
[KDeta,KDtheta,KReta,KRtheta] = stiffmats(EImat,Lmat);
stiffmat = [KDeta KDtheta;KReta KRtheta];

%----------------------
% Make the global stiffness matrix (static reduction of elements with zero
% inertia)
Kmat = KDtheta*inv(KRtheta)*(KReta);
Kmat = KDeta - Kmat;
%Kmat = [KDeta KDtheta;KReta KRtheta];

%----------------------
% Define the first-order term
Dmat = zeros(size(Kmat));

%-----------------------
% Reduce to a set of first order differential equations
M = diag(Mmat);
invM = diag(1./Mmat);
%invM = diag([1./Mmat;ones(size(Mmat))]);
%H = [zeros(size(Kmat)) eye(size(Kmat));invM*Kmat invM*Dmat];
[Veigs,Deigs] = eigs(Kmat,M,6,'LR');
Deigs = diag(Deigs);

nels = length(EImat);
VDisp = Veigs([1:nels],:);

VVelo = Veigs([nels+1:2*nels]);
