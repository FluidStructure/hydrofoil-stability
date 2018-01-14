function [massmat,damping,stiffness] = addHydroIntegralsVertical(massmat,damping,stiffness,F,i)

i = ((i-1)*6) + 1;

% Added mass in the transverse (x-direction)
massmat(i+0,i+0) = massmat(i+0,i+0) + F(2,4);
massmat(i+0,i+1) = massmat(i+0,i+1) + F(2,1);
massmat(i+0,i+5) = massmat(i+0,i+5) - F(2,7);
% Added mass in the for-aft (y-direction)
massmat(i+1,i+0) = massmat(i+1,i+0) + F(1,4);
massmat(i+1,i+1) = massmat(i+1,i+1) + F(1,1);
massmat(i+1,i+5) = massmat(i+1,i+5) + F(1,7);
% Added mass in the rotational direction
massmat(i+5,i+0) = massmat(i+5,i+0) - F(3,4);
massmat(i+5,i+1) = massmat(i+5,i+1) + F(3,1);
massmat(i+5,i+5) = massmat(i+5,i+5) + F(3,7);

% Damping in the transverse (x-direction)
damping(i+0,i+0) = damping(i+0,i+0) + F(2,5);
damping(i+0,i+1) = damping(i+0,i+1) + F(2,2);
damping(i+0,i+5) = damping(i+0,i+5) - F(2,8);
% Damping in the for-aft (y-direction)
damping(i+1,i+0) = damping(i+1,i+0) + F(1,5);
damping(i+1,i+1) = damping(i+1,i+1) + F(1,2);
damping(i+1,i+5) = damping(i+1,i+5) + F(1,8);
% Damping in the rotational direction
damping(i+5,i+0) = damping(i+5,i+0) - F(3,5);
damping(i+5,i+1) = damping(i+5,i+1) + F(3,2);
damping(i+5,i+5) = damping(i+5,i+5) + F(3,8);

% Stiffness in the transverse (x-direction)
stiffness(i+0,i+0) = stiffness(i+0,i+0) + F(2,6);
stiffness(i+0,i+1) = stiffness(i+0,i+1) + F(2,3);
stiffness(i+0,i+5) = stiffness(i+0,i+5) - F(2,9);
% Stiffness in the for-aft (y-direction)
stiffness(i+1,i+0) = stiffness(i+1,i+0) + F(1,6);
stiffness(i+1,i+1) = stiffness(i+1,i+1) + F(1,3);
stiffness(i+1,i+5) = stiffness(i+1,i+5) + F(1,9);
% Stiffness in the rotational direction
stiffness(i+5,i+0) = stiffness(i+5,i+0) - F(3,6);
stiffness(i+5,i+1) = stiffness(i+5,i+1) + F(3,3);
stiffness(i+5,i+5) = stiffness(i+5,i+5) + F(3,9);