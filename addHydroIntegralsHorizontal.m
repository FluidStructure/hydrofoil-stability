function [massmat,damping,stiffness] = addHydroIntegralsHorizontal(massmat,damping,stiffness,F,i)

i = ((i-1)*6) + 1;

% Added mass in the vertical (z-direction)
massmat(i+2,i+2) = massmat(i+2,i+2) + F(2,4);
massmat(i+2,i+1) = massmat(i+2,i+1) + F(2,1);
massmat(i+2,i+3) = massmat(i+2,i+3) + F(2,7);
% Added mass in the for-aft (y-direction)
massmat(i+1,i+2) = massmat(i+1,i+2) + F(1,4);
massmat(i+1,i+1) = massmat(i+1,i+1) + F(1,1);
massmat(i+1,i+3) = massmat(i+1,i+3) + F(1,7);
% Added mass in the rotational direction
massmat(i+3,i+2) = massmat(i+3,i+2) + F(3,4);
massmat(i+3,i+1) = massmat(i+3,i+1) + F(3,1);
massmat(i+3,i+3) = massmat(i+3,i+3) + F(3,7);

% Damping in the transverse (z-direction)
damping(i+2,i+2) = damping(i+2,i+2) + F(2,5);
damping(i+2,i+1) = damping(i+2,i+1) + F(2,2);
damping(i+2,i+3) = damping(i+2,i+3) + F(2,8);
% Damping in the for-aft (y-direction)
damping(i+1,i+2) = damping(i+1,i+2) + F(1,5);
damping(i+1,i+1) = damping(i+1,i+1) + F(1,2);
damping(i+1,i+3) = damping(i+1,i+3) + F(1,8);
% Damping in the rotational direction
damping(i+3,i+2) = damping(i+3,i+2) + F(3,5);
damping(i+3,i+1) = damping(i+3,i+1) + F(3,2);
damping(i+3,i+3) = damping(i+3,i+3) + F(3,8);

% Stiffness in the transverse (z-direction)
stiffness(i+2,i+2) = stiffness(i+2,i+2) + F(2,6);
stiffness(i+2,i+1) = stiffness(i+2,i+1) + F(2,3);
stiffness(i+2,i+3) = stiffness(i+2,i+3) + F(2,9);
% Stiffness in the for-aft (y-direction)
stiffness(i+1,i+2) = stiffness(i+1,i+2) + F(1,6);
stiffness(i+1,i+1) = stiffness(i+1,i+1) + F(1,3);
stiffness(i+1,i+3) = stiffness(i+1,i+3) + F(1,9);
% Stiffness in the rotational direction
stiffness(i+3,i+2) = stiffness(i+3,i+2) + F(3,6);
stiffness(i+3,i+1) = stiffness(i+3,i+1) + F(3,3);
stiffness(i+3,i+3) = stiffness(i+3,i+3) + F(3,9);