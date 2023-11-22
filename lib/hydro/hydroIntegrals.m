function [F] = hydroIntegrals(xc,yc,CL,NACA,HALF_NPANELS,rho,Uinf,Ls)

%
% A function that computes the hydrodynamic integrals for displacment and
% rotation of a foil.
%
% Inputs:
%   xc,yc = centre-point about which to calculate hydrodynamic-moment integrals
%   NACA = 4-digit NACA number
%   rho = density of the fluid
%
% Outputs:
%
%   The matrix F which can be used to give Fx, Fy, and Fr when:
%       [Fx;Fy;Fr] = F*[x_ddot;x_dot;x;y_ddot;y_dot;y;r_ddot;r_dot;r]
%
%   Xi,Xd,Xs = Hydrodynamic integrals of forcing in X(horizontal) direction for 
%   inertia (eta_ddot), damping (eta_dot) and stiffness (eta) respectively
%
%   Yi,Yd,Ys = Hydrodynamic integrals of forcing in Y(vertical) direction for 
%   inertia (eta_ddot), damping (eta_dot) and stiffness (eta) respectively
%
%   Ri,Rd,Rs = Hydrodynamic integrals of moments for inertia (theta_ddot), damping
%   (theta_dot) and stiffness (theta) respectively
%

%--------------------------
% GENERATE THE AIRFOIL SHAPE

% C_L = CL;                               % Chord length
% HALF_NPANELS = HALF_PANELS;             % No. of panels on one side of aerofoil
% MM = floor(mod(NACA/1000,10))*CL/100;   % Camber (first digit of 4-digit NACA)
% P = floor(mod(NACA/100,10))*CL/100;     % Position of centre of pressure (% of chord)
% T = floor(mod(NACA/1,100))*CL/100;      % Thickness

% Get NACA section using function aerofoil_geom
% [PANEND,PLENGT,UVNORM,UVTANG,CP] = aerofoil_geom(HALF_NPANELS,C_L,MM,P,T);

% Get NACA section using function NACAsec
% [PANEND,PLENGT,CP] = NACAsec(HALF_NPANELS,HALF_NPANELS,NACA,4);
% PANEND = PANEND.*CL;
% PLENGT = PLENGT.*CL;
% PLENGT=PLENGT';

% Get NACA section using function naca4gen
iaf.designation=num2str(NACA);
iaf.n=HALF_NPANELS;
iaf.HalfCosineSpacing=1;
iaf.wantFile=0;
iaf.datFilePath='./'; % Current folder
iaf.is_finiteTE=0;
af = naca4gen(iaf);
PANEND(:,1) = flipud(af.x).*CL;
PANEND(:,2) = flipud(af.z).*CL;
CP(:,1) = (PANEND([1:(length(PANEND)-1)],1) + PANEND([2:(length(PANEND))],1))./2;
CP(:,2) = (PANEND([1:(length(PANEND)-1)],2) + PANEND([2:(length(PANEND))],2))./2;

%--------------------------
% GENERATE ELEMAT AND THE INFLUENCE MATRICES
elemat(:,1) = CP(:,1);
elemat(:,2) = CP(:,2);
nn = size(PANEND,1);
dy = PANEND([2:nn],2) - PANEND([1:nn-1],2);
dx = PANEND([2:nn],1) - PANEND([1:nn-1],1);
elemat(:,3) = atan2(dy,dx);
elemat(:,4) = sqrt(dx.^2 + dy.^2);
elemat(:,5) = zeros(size(elemat(:,1)));

[ua,ub,va,vb] = pvt_lin_finic_xy(elemat);

%--------------------------
% Get the normal and tangential influence coefficients between panels
for i = 1:size(ua,1)
    Na(i,:) = va(i,:).*cos(elemat(i,3)) - ua(i,:).*sin(elemat(i,3));
    Nb(i,:) = vb(i,:).*cos(elemat(i,3)) - ub(i,:).*sin(elemat(i,3));
    
    Ta(i,:) = va(i,:).*sin(elemat(i,3)) + ua(i,:).*cos(elemat(i,3));
    Tb(i,:) = vb(i,:).*sin(elemat(i,3)) + ub(i,:).*cos(elemat(i,3));
end

%Na = Na*-1;Nb = Nb*-1;

% Normal influence coefficient
ne = size(Na,1);
N = zeros(ne+1,ne+1);
N([1:ne],[1:ne]) = N([1:ne],[1:ne]) + Na;
N([1:ne],[2:ne+1]) = N([1:ne],[2:ne+1]) + Nb;

% Tangential influence coefficent
T = zeros(ne,ne+1);
T([1:ne],[1:ne]) = T([1:ne],[1:ne]) + Ta;
T([1:ne],[2:ne+1]) = T([1:ne],[2:ne+1]) + Tb;


% BOUNDARY CONDITIONS TO NORMAL INFLUENCE MATRIX!!!
% Apply continuity boundary condition at first-last panel
%N(ne+1,1)=1;N(ne+1,ne+1)=-1;
% Apply Kutta condition at trailing edge
%N(ne+1,1) = 1;
N(ne+1,1) = 1;N(ne+1,ne+1)=1;
% Apply total bound vorticity equals zero
%N(ne+1,:) = ones(size(N(ne+1,:)));
% Apply total bound vorticity equals free-stream vorticity
% N(ne+1,:) = zeros(size(N(ne+1,:)));
% for i=1:ne
% 	N(ne+1,i) = N(ne+1,i) + 0.5*PLENGT(i);
% 	N(ne+1,i+1) = N(ne+1,i+1) + 0.5*PLENGT(i);
% end

% Generate Velocity Potential Influence coefficient
% (integration of tangential velocity)
PHI = zeros(size(T));
phi_node1(1,:) = zeros(1,size(T,2));
for i = 1:(ne)
    %phi_node2 = phi_node1 - (Uinf.*cos(elemat(i,3)) + T(i,:)).*elemat(i,4);
    phi_node2 = phi_node1 - T(i,:).*elemat(i,4);
    PHI(i,:) = (phi_node1 + phi_node2)./2;
    %PHI(i,:) = 0*(phi_node1 + phi_node2)./2;   % Debugging only
    
    phi_node1 = phi_node2;
end

% Get the inverse of the normal velocity matrix
invN = inv(N);
Lx = elemat(:,1) - xc;
Ly = elemat(:,2) - yc;
ds = elemat(:,4);

% Evaluate the hydrodynamic integrals for PRESSURE
K = (Uinf.*cos(elemat(:,3))) + (T*invN)*([Uinf.*sin(elemat(:,3));0]);

P(:,1) = (rho)*(PHI*invN)*[sin(elemat(:,3));0];
P(:,2) = (-1*rho*K).*((T*invN)*([sin(elemat(:,3));0]) + cos(elemat(:,3)));
P(:,3) = zeros(size(P(:,1)));

P(:,4) = (rho)*(PHI*invN)*[cos(elemat(:,3));0];
P(:,5) = (-1*rho*K).*((T*invN)*([cos(elemat(:,3));0]) + sin(elemat(:,3)));
%P(:,5) = P(:,5)*-1;
P(:,6) = zeros(size(P(:,1)));

P(:,7) = (rho)*(PHI*invN)*([Lx.*cos(elemat(:,3));0]) + (rho)*(PHI*invN)*([Ly.*sin(elemat(:,3));0]);
P(:,8) = (-1*rho*K).*(T*invN*[Lx.*cos(elemat(:,3));0] + Lx.*sin(elemat(:,3))) + (-1*rho*K).*((T*invN*[Ly.*sin(elemat(:,3));0]) - Ly.*cos(elemat(:,3))) + (rho)*(PHI*invN)*([Uinf.*cos(elemat(:,3));0]);
P(:,9) = (-1*rho*K).*(T*invN*[Uinf.*cos(elemat(:,3));0]) + (rho*K).*(Uinf.*sin(elemat(:,3)));

P(:,10) = (0.5*rho).*(Uinf.^2 - K.*K);       % Constant Pressure Term

% Evaluate the hydrodynamic integrals for forcing in various directions
% x-direction
TMP = zeros(size(P));
for i = 1:size(P,2)
    TMP(:,i) = (Ls.*ds.*sin(elemat(:,3))).*P(:,i).*0;
end
F(1,:) = sum(TMP,1);

% y-direction
TMP = zeros(size(P));
for i = 1:size(P,2)
    TMP(:,i) = (-1*Ls.*ds.*cos(elemat(:,3))).*P(:,i);
end
F(2,:) = sum(TMP,1);

% rotational-direction
TMP = zeros(size(P));
for i = 1:size(P,2)
    TMP(:,i) = (-1*Ls.*ds.*cos(elemat(:,3))).*Lx.*P(:,i) + (-1*Ls.*ds.*sin(elemat(:,3))).*Ly.*P(:,i);
end
F(3,:) = sum(TMP,1);

%figure;plot(PANEND(:,1),PANEND(:,2),'rx-')
%figure;gam = invN*([Uinf.*sin(elemat(:,3));0]);plot(gam,'bx-');grid

