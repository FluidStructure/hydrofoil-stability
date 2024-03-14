%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [G] = getfr(In,Idim,On,Odim,H,Fmin,Fmax,Fn,numberNodes,invindmat)
%[G] = getfr(Nv-1,1,Nv-1,1,H,0.1,1000,500,numberNodes,invindmat);

% Function that gets the frequency response from a first order system with
% state matrix H

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lh = size(H,1); B = zeros(lh,1); C = zeros(1,lh); D = 0;

% Get a reference index for nodepositions and dimensions within H
Hrows = [1:(size(H,1)/2)]';
Hind = zeros((numberNodes-1)*6,1); Hind(invindmat) = Hrows;
Inn = Hind((In-1)*6 + Idim); Onn = Hind((On-1)*6 + Odim);

if (Inn == 0)||(Onn == 0)
    error('Error looking up node and dimension within H')
end

A = H; B(Inn,1) = 1; C(1,Onn) = 1;
F = logspace(log10(Fmin),log10(Fmax),Fn);
W = (2*pi*F)*1i;
keyboard
[G] = ltifreq(A,B,C,D,W);

% ey = eye (size (A)); lw = length (W); G = zeros (lw, 1);
% for ii = 1:lw
%     out = (W(ii)*ey-A)\B;
%     out = C*out + D;
%     G(ii) = out;
% end

%K = H([(lh/2)+1:lh],[1:(lh/2)]);
%E = zeros(size(K,1),1);P = E';
%E(Inn) = 1;P(Onn) = 1;
%dst = P*(K\E);
dst = abs(G(1)); % Get the static deflection factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot a bode diagram

figure;
subplot(2,1,1);
loglog(F,(abs(G)./dst),'k')
grid; ylabel('Amplitude ratio')

subplot(2,1,2);
semilogx(F,angle(G),'k')
grid; xlabel('Forcing Frequency (Hz)'); ylabel('Phase shift')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
