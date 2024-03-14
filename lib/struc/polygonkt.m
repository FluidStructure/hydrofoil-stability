%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ixx,Iyy,A,Properties of_airfoil section%%%%%%
function [Ixx,Iyy,xc,yc,A]=polygonkt(xnaca,ynaca,chord)
N_points=size(xnaca,1);
znaca_orig(1:N_points,1)= xnaca + 1i*ynaca; %(m)

 A=0;
 sumx=0;
 sumy=0;
for ip=2:N_points
dx=abs(    ( real(znaca_orig(ip)) + chord/2 ) - ( real( znaca_orig(ip-1)) + chord/2 )     ); 

dA= ( abs(imag(znaca_orig(ip))) + abs(imag( znaca_orig(ip-1))) )  *  dx / 2;

xp = ( real(znaca_orig(ip)) + chord/2  + real(znaca_orig(ip-1)) + chord/2 ) / 2; 

yp = ( imag(znaca_orig(ip))/2 + imag(znaca_orig(ip-1))/2 ) / 2;

sumx = sumx + dA*xp;

sumy = sumy + dA*yp;

A=A+dA;

end
xc=sumx/A;
yc=sumy/A;
xc_chord=xc/chord;
yc_chord=yc/chord;
A_chord=A/chord^2;


sumx=0;
for ip=2:N_points
dx=abs(    ( real(znaca_orig(ip)) + chord/2 ) - ( real( znaca_orig(ip-1)) + chord/2 )     );

dA= ( abs(imag(znaca_orig(ip))) + abs(imag( znaca_orig(ip-1))) )  *  dx / 2;

xps = ( ( real(znaca_orig(ip)) + chord/2  - xc  + real(znaca_orig(ip-1)) + chord/2  - xc ) / 2 )^2;

sumx = sumx + dA*xps;
end

Iyy=sumx/chord^4;




sumy=0;
A=0;
for ip=2:N_points
dy = abs( imag(znaca_orig(ip)) - imag(znaca_orig(ip-1)) ); 

dA=( abs( real(znaca_orig(ip)) + chord/2 - xc ) + abs( real(znaca_orig(ip-1))+ chord/2 - xc ) ) * dy / 2;

yps = ((imag(znaca_orig(ip))+imag(znaca_orig(ip-1)))/2  - yc)^2;

sumy = sumy + dA*yps;

A=A+dA;

end

Ixx=sumy/chord^4;
%A_chord=A/chord^2



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ixx,Iyy,A,Properties of_airfoil section%%%%%%




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ixx,Iyy,A,Properties of_airfoil section%%%%%%
% 
%  A=0;
%  sumx=0;
%  sumy=0;
% for ip=2:N_points
% dx=abs(    (real(zet(ip)) + chord/2 ) - ( real(zet(ip-1)) + chord/2 )     ); 
% 
% dA= ( abs(imag(zet(ip))) + abs(imag(zet(ip-1))) )  *  dx / 2;
% 
% xp = ( real(zet(ip)) + chord/2  + real(zet(ip-1)) + chord/2 ) / 2; 
% 
% yp = ( imag(zet(ip)) + imag(zet(ip-1)) ) / 2;
% 
% sumx = sumx + dA*xp;
% 
% sumy = sumy + dA*yp;
% 
% A=A+dA;
% 
% end
% xc=sumx/A;
% yc=sumy/A;
% xc_chord=xc/chord
% yc_chord=yc/chord
% A_chord=A/chord^2
% 
% 
% 
% sumx=0;
% for ip=2:N_points
% dx=abs(    ( real(zet(ip)) + chord/2 ) - ( real(zet(ip-1)) + chord/2 )     );
% 
% dA= ( abs(imag(zet(ip))) + abs(imag(zet(ip-1))) )  *  dx / 2;
% 
% xps = ( ( real(zet(ip)) + chord/2  - xc  + real(zet(ip-1)) + chord/2  - xc ) / 2 )^2;
% 
% sumx = sumx + dA*xps;
% end
% 
% Iyy=sumx/chord^4
% 
% 
% 
% sumy=0;
% A=0;
% for ip=2:N_points
% dy = abs( imag(zet(ip)) - imag(zet(ip-1)) ); 
% 
% dA=( abs( real(zet(ip)) + chord/2 - xc ) + abs( real(zet(ip-1))+ chord/2 - xc ) ) * dy / 2;
% 
% yps = (imag(zet(ip)) - yc)^2;
% 
% sumy = sumy + dA*yps;
% 
% A=A+dA;
% 
% end
% 
% Ixx=sumy/chord^4
% %A_chord=A/chord^2
% 
% 
% 
% 
% 
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ixx,Iyy,A,Properties of_airfoil section%%%%%%