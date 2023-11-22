% stvt_04
%
% This is the starting vortex problem solved using linear vortex elements
% around the aerofoil and shedding discrete vortex blobs that are equal in
% strength to the change in bound vorticity on the surface of the aerofoil
% at each timestep.
%
% This version 4 of the code does not shed blobs and just performs a static
% potential flow calculation.  This is so that the blob-shedding codes have
% a reference.
%

%close all;clear all;

%--------------------------
% Some calculation parameters
Uinfm = 5;  % Mean flow velocity (m/s)
alph = 5*pi/180;  % Angle of attack (radians)

ntstp = 300; % Number of time-steps
dt = 10e-3;

csize = 0.01;

FUDGE = 1;
%--------------------------
% GENERATE THE AIRFOIL SHAPE

C_L = 1.0 % Chord length
HALF_NPANELS = 60 %No. of panels on one side of aerofoil
MM = 0.02 % Camber
P = 0.4 % Position of centre of pressure (% of chord)
T = 0.12 % Thickness

% Get NACA section using function aerofoil_geom
% [PANEND,PLENGT,UVNORM,UVTANG,CP] = aerofoil_geom(HALF_NPANELS,C_L,MM,P,T);

% Get NACA section using function NACAsec
[PANEND,PLENGT,CP] = NACAsec(HALF_NPANELS,HALF_NPANELS,2412,4);
PLENGT=PLENGT';

% PANEND - panel end points
% PLENGT - panel lengths
% UVNORM UVTANG - unit vectors (i,j)
% CP - control point co-ordinates (x,y)

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

% %--------------------------
% % FOR DEBUGGING ONLY
% % GET THE INFLUENCE OF EACH PANEL ON ITSELF
% ne = size(elemat,1);
% for i = 1:ne
%     Us(i) = ua(i,i) + ub(i,i);
%     Vs(i) = va(i,i) + vb(i,i);
% end
% 
% %return

%--------------------------
% Get the normal and tangential influence coefficients between panels
for i = 1:size(ua,1)
    Na(i,:) = -va(i,:).*cos(elemat(i,3)) + ua(i,:).*sin(elemat(i,3));
    Nb(i,:) = -vb(i,:).*cos(elemat(i,3)) + ub(i,:).*sin(elemat(i,3));
    
    Ta(i,:) = va(i,:).*sin(elemat(i,3)) + ua(i,:).*cos(elemat(i,3));
    Tb(i,:) = vb(i,:).*sin(elemat(i,3)) + ub(i,:).*cos(elemat(i,3));
end

ne = size(Na,1);
N = zeros(ne+1,ne+1);
N([1:ne],[1:ne]) = N([1:ne],[1:ne]) + Na;
N([1:ne],[2:ne+1]) = N([1:ne],[2:ne+1]) + Nb;
% Apply continuity boundary condition at first-last panel
%N(ne+1,1)=1;N(ne+1,ne+1)=-1;
% Apply Kutta condition at trailing edge
N(ne+1,HALF_NPANELS+1) = 1;
% Apply total bound vorticity equals zero
%N(ne+1,:) = ones(size(N(ne+1,:)));


invN = inv(N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% BEGIN TIME-STEPPING
totboundvort = 0;
vormat = [];
%fig1 = figure;
xax = C_L+(Uinfm*dt*ntstp)*1.1;
cntr = 1;
ctime = dt;
Uinfmat = zeros([1:ntstp],1);
for ts = 1:ntstp
    
    % Ramp-up Uinf over the first N seconds
    Nts = ntstp*dt;		% Set to -1 to impulsively start the flow (straight up to Uinfm)
	% Nts = -1;
    if ctime < Nts;
        % Uinf = Uinfm*(ctime/Nts);			% Ramp-up the flow speed linearly
		Uinf = Uinfm*(1-exp(-2*ctime));		% Ramp-up flow speed exponentially
    else
        Uinf = Uinfm;
    end
    Uinfmat(ts,1) = Uinf;	

    % Get the influence of the shed vorticies on the aerofoil
    [Nsv] = vort_on_foil(vormat,elemat);
    % solve for the vortex element strengths
    RHS = (Uinf*cos(alph)).*sin(elemat(:,3)) - (Uinf*sin(alph)).*cos(elemat(:,3));
    RHSsv = RHS + Nsv;RHSbc = -1*RHSsv;RHSbc(ne+1,1) = 0;gam = invN*RHSbc;
    
    % Get the change in total bound vorticity
    nn = size(gam,1);
    gamav = (gam(1:nn-1) + gam(2:nn))./2;
    tmp1 = sum(gamav.*PLENGT');
    Dgam = tmp1 - totboundvort;
    totboundvort = tmp1;
    
    %--------------------------
    % Create some graphical output

    %figure(fig1)
if cntr == 1
    fig1 = figure;
    plot(PANEND(:,1),PANEND(:,2),'k-x')
    %axis([0 xax -xax/2 xax/2])
    axis equal    
    hold on;
    if isempty(vormat)==0
        plot(vormat(:,1),vormat(:,2),'r.')
    end
    plot(elemat(:,1),elemat(:,2),'ko')
    hold off;
    grid
    title('NACA 2412 section profile (x=node points, o=collocation point)')
    xlabel('x-coordinate (m)')
    ylabel('y-coordinate (m)')
    pause(0.001)

%   print(fig1,'-deps','NACAsec')
end

disp(['Time step = ' num2str(ts)])

     Vt = Ta*gam(1:ne) + Tb*gam(2:ne+1) + (Uinf*cos(alph)).*cos(elemat(:,3)) + (Uinf*sin(alph)).*sin(elemat(:,3));
     Vn = Na*gam(1:ne) + Nb*gam(2:ne+1) + (Uinf*cos(alph)).*sin(elemat(:,3)) - (Uinf*sin(alph)).*cos(elemat(:,3));
     
     U = ua*gam(1:ne) + ub*gam(2:ne+1) + Uinf*cos(alph);
     V = va*gam(1:ne) + vb*gam(2:ne+1) + Uinf*sin(alph);
     
     %fig1 = figure;
     %quiver(elemat(:,1),elemat(:,2),U,V);
     %axis equal
     %grid


	Cp = 1 - (Vt/Uinf).^2;	% Coefficient of Pressure
	xnd = elemat(:,1)/C_L;		% Non-dimensional x-coordinate
%	fig2 = figure;
%	plot(xnd,Cp,'b-o');

% Get the Unsteady Part of the Bernoulli equation



	% Calculate the coefficient of lift (and drag)
	Cf = Cp.*PLENGT';
	Cl = -1*sum(Cf.*cos(elemat(:,3)));

	Cpmat(:,cntr) = [Cp];
	Clmat(cntr,:) = [ctime Cl];

	% save 10deg-mark.mat xnd Cp
	% return

	cntr = cntr + 1;
    ctime = ctime + dt;
      
end

