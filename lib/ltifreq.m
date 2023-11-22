function [G] = ltifreq(A,B,C,D,W)

%
% Calculates the frequency response from an input in state space form
%

ey = eye (size (A));
lw = length (W);
G = zeros (lw, 1);

for ii = 1:lw,
    out = (W(ii)*ey-A)\B;
    out = C*out + D;
    G(ii) = out;
end

