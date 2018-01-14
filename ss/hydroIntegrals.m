function [mt,dt,kt,mr,dr,kr] = hydroIntegrals(alpha,C_L,HALF_NPANELS,MM,P,T)

%
% This script calculates the mean hydrodynamic integrals of mass, damping and stiffness for 
% translational and rotational motion (about some defined centre point) for the NACA airfoil
% with given properties.
%
% NOTE:
% The NACA four-digit wing sections define the profile by:[1]
%   1. One digit describing maximum camber as percentage of the chord.
%   2. One digit describing the distance of maximum camber from the airfoil leading edge in tens of percents of the chord.
%   3. Two digits describing maximum thickness of the airfoil as percent of the chord.
%
% Inputs:
% alpha = mean angle of attack of airfoil (radians)
% C_L = 1.0 % Chord length
% HALF_NPANELS = 60 %No. of panels on one side of aerofoil
% MM = 0.02 % Camber
% P = 0.4 % Position of centre of pressure (% of chord)
% T = 0.12 % Thickness
% NOTE: Example above is for a NACA 2412 profile 
%
% Outputs:
% mt,dt,kt = hydrodynamic integrals of added mass, damping and stiffness respectively for translational motion
% mr,dr,kr = dydrodynamic integrals of added mass, damping and stiffness respectively for rotational motion
%

%--------------------------
% Get NACA section using function NACAsec
[PANEND,PLENGT,CP] = NACAsec(HALF_NPANELS,HALF_NPANELS,2412,4);
PLENGT=PLENGT';
% PANEND - panel end points
% PLENGT - panel lengths
% CP - control point co-ordinates (x,y)

%--------------------------
% GENERATE ELEMAT
elemat(:,1) = CP(:,1);
elemat(:,2) = CP(:,2);
nn = size(PANEND,1);
dy = PANEND([2:nn],2) - PANEND([1:nn-1],2);
dx = PANEND([2:nn],1) - PANEND([1:nn-1],1);
elemat(:,3) = atan2(dy,dx);
elemat(:,4) = sqrt(dx.^2 + dy.^2);
elemat(:,5) = zeros(size(elemat(:,1)));

%--------------------------
% GENERATE THE INFLUENCE MATRICES

% Velocity Potential (integrate around the airfoil from the trailing edge)
[ua,ub,va,vb] = pvt_lin_finic_xy(elemat);







