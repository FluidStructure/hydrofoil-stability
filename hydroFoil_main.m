%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Hydroelastic Stability Analysis of an Inverted T Hydrofoil
% M. Pitman 2009
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialise

clear; close all

if ~exist('getfr.m','file')
    addpath(genpath('lib')); addpath(genpath('data'));
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data

dataflg=0;

switch dataflg
    case{0}; data_orig;
    case{1}; data_lelong;
end

% Number of iterations and eigenmodes to return & the eigenmode to plot
numits = 300; numeigs = 20; NMplot = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% main code

cntr = 1;
for Uinf = Uinfm
    disp(['Uinf = ' num2str(Uinf)]);

    stiffness_main % 3D stiffness matrix for an inverted T-shaped hydrofoil

    hydro_main % Generate the hydrodynamic integrals for each node

    fsi_main % couple fluid and structure & solve system

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
if (length(Uinfm)==1);  plot_bode; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
