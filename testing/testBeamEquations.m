%
% Little script that plots the SF, BM, Slope and Displacment curves for 
% a simply supported beam under some boundary conditions
%

P1 = 1;             % P1 (for clamped-clamped beam with zero gradient at each end)
P2 = 1;             % P2 (for clamped-clamped beam with some gradient at RHS end)

EI = 0.1;           % Properties of the beam

L = 1;              % Length of the beam
x = [0:L/100:L];    % x-coordinates to plot the curves over

%---------------------
% Case 1 - clamped-clamped beam with zero gradient at each end
P = P1;
SF1 = P.*ones(size(x));
BM1 = P.*x - (P*L/2);
theta1 = (1/EI).*((P/2).*(x.^2) - (P*L/2).*x);
delta1 = (1/EI).*((P/6).*(x.^3) - (P*L/4).*(x.^2));

figure;hold on;
plot(x,[SF1;BM1;theta1;delta1])
hold off
legend('SF','BM','theta','delta','Location','NorthWest')
grid
title('Clamped-clamped beam with zero gradient at each end')

%---------------------
% Case 2 - clamped-clamped beam with some gradient at RHS end
P = P2;
SF1 = P.*ones(size(x));
BM1 = P.*x - (P*L/3);
theta1 = (1/EI).*((P/2).*(x.^2) - (P*L/3).*x);
delta1 = (1/EI).*((P/6).*(x.^3) - (P*L/6).*(x.^2));

figure;hold on;
plot(x,[SF1;BM1;theta1;delta1])
hold off
legend('SF','BM','theta','delta','Location','NorthWest')
grid
title('Clamped-clamped beam with zero gradient at each end')


%--------------------
% Superpose results given displacment and angle at RHS
delta = -0.5;  % Unit displacement
theta = 1;  % 45degree angle

Pdisp = -1*(delta.*EI)*12./(L.^3);
Pangl = (theta.*EI)*6./(L.^2);

P = Pdisp;
SF1 = P.*ones(size(x));
BM1 = P.*x - (P*L/2);
theta1 = (1/EI).*((P/2).*(x.^2) - (P*L/2).*x);
delta1 = (1/EI).*((P/6).*(x.^3) - (P*L/4).*(x.^2));

P = Pangl;
SF2 = P.*ones(size(x));
BM2 = P.*x - (P*L/3);
theta2 = (1/EI).*((P/2).*(x.^2) - (P*L/3).*x);
delta2 = (1/EI).*((P/6).*(x.^3) - (P*L/6).*(x.^2));

SF = SF1 + SF2;
BM = BM1 + BM2;
theta = theta1 + theta2;
delta = delta1 + delta2;

figure;hold on;
plot(x,[SF;BM;theta;delta])
hold off
legend('SF','BM','theta','delta','Location','NorthWest')
grid
title('Clamped-clamped beam with displacement and rotation at RHS')
axis equal

Prhs = Pdisp + Pangl;
Mrhs = (Pdisp.*L/2) + (2.*Pangl.*L/3);
disp(['Force at RHS = ' num2str(Prhs)])
disp(['Moment at RHS = ' num2str(Mrhs)])
