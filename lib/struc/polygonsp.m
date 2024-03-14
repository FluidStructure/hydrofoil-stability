function sp = polygonsp(x,y)
%
% Computes the various section properties (area moments) of a planar polygon 
% that is defined by its vertices. The polygon can be simple, or it can consist
% of disconnected multiple polygons. The area to be integrated is defined by the
% sequence of vertices contained in the vectors x and y. The user must ensure 
% that first and last vertices are the same.
%
% The function computes the following section properties and returns them in
% the structure sp:
%
% sp.Area              = area
% sp.MAx               = first moment of area about the x axis
% sp.MAy               = first moment of area about the y axis
% sp.Ixx               = second moment of area about the x axis
% sp.Iyy               = second moment of area about the y axis
% sp.Izz               = second moment of area about the z axis
% sp.Ixy               = cross product of area about the x and y axes
% sp.Xrog              = radius of gyration about the x axis
% sp.Yrog              = radius of gyration about the y axis
% sp.Zrog              = radius of gyration about the z axis
% sp.Xc                = x location of the centroid of the cross section
% sp.Yc                = x location of the centroid of the cross section
% sp.Icxx              = second moment of area about the axis parallel to the x axis 
%                        and passing through the centroid
% sp.Icyy              = second moment of area about the axis parallel to the y axis 
%                        and passing through the centroid
% sp.Iczz              = second moment of area about the axis parallel to the x axis 
%                        and passing through the centroid
% sp.Icxy              = cross product of area about the axes parallel to the x and
%                        y axes and passing through the centroid
% sp.Xcrog             = radius of gyration about the axis parallel to the x axis
%                        and passing through the centroid
% sp.Ycrog             = radius of gyration about the axis parallel to the y axis
%                        and passing through the centroid
% sp.Zcrog             = radius of gyration about the axis parallel to the z axis
%                        and passing through the centroid
% sp.ImaxPrincipalAxes = maximum principal second moment of area
% sp.IminPrincipalAxes = minimum principal second moment of area
% sp.MajorAxisRotDeg   = angle of rotation (degrees) of the major principal axis
%                        relative to the x axis. Positive in an anticlockwise direction
%
% References:
%
% Marin, J. (1984). Computing columns, footings and gates through moments
% of area. Computers & Structures, Vol. 18, No. 2, pp. 343–349.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Check the number of input arguments and proceed accordingly.
if nargin == 0
  runTestCases();
  return
elseif nargin ~= 2
  error('Incorrect number of arguments supplied to polysectprops')
end

P=size(x,1)-1; i   = 1:P; ip1 = i+1;

xi   = x(i);xip1 = x(ip1);yi   = y(i);yip1 = y(ip1);

wi = xi.*yip1 - xip1.*yi; % Equation (2.4) from Marin.

% Eqns (2.7) to (2.12) from Marin and other known results.
sp.Area  = sum(wi)/2;
sp.MAx   = sum(wi.*(yi + yip1))/6;
sp.MAy   = sum(wi.*(xi + xip1))/6;
sp.Ixx   = sum(wi.*((yi + yip1).^2 - yi.*yip1))/12;
sp.Iyy   = sum(wi.*((xi + xip1).^2 - xi.*xip1))/12;
sp.Izz   = sp.Ixx + sp.Iyy;
sp.Ixy   = sum(wi.*((xi + xip1).*(yi + yip1) + xi.*yi + xip1.*yip1))/24;
sp.Xrog  = sqrt(sp.Ixx/sp.Area);
sp.Yrog  = sqrt(sp.Iyy/sp.Area);
sp.Zrog  = sqrt(sp.Izz/sp.Area);
sp.Xc    = sp.MAy/sp.Area;
sp.Yc    = sp.MAx/sp.Area;
sp.Icxx  = sp.Ixx - sp.Yc^2*sp.Area;
sp.Icyy  = sp.Iyy - sp.Xc^2*sp.Area;
sp.Iczz  = sp.Icxx + sp.Icyy;
sp.Icxy  = sp.Ixy - sp.Xc*sp.Yc*sp.Area;
sp.Xcrog = sqrt(sp.Icxx/sp.Area);
sp.Ycrog = sqrt(sp.Icyy/sp.Area);
sp.Zcrog = sqrt(sp.Iczz/sp.Area);

Iavg = (sp.Icxx + sp.Icyy)/2;
Idif = (sp.Icxx - sp.Icyy)/2;
sp.ImaxPrincipalAxes = Iavg + sqrt(Idif^2 + sp.Icxy^2);
sp.IminPrincipalAxes = Iavg - sqrt(Idif^2 + sp.Icxy^2);

% Compute the angle of the major axis of the section relative to the x axis.
% Here the major (strong) and minor (weka) axes pass through the centroid of
% the section. Handle cases including those where the section is symmetrical
% about the x and y axes, which results in Icxy = 0 and Ixx - Iyy = 0.
maxIcxxIcyy = max(sp.Icxx,sp.Icyy);
num = 2*sp.Icxy/maxIcxxIcyy; den = (sp.Icyy-sp.Icxx)/maxIcxxIcyy;
eps = 1.0e-08;
if abs(num) <= eps;  num = 0; end
if abs(den) <= eps;  den = 0; end
thetad = atan2d(num,den)/2;
Iee = sp.Icxx*cosd(thetad)^2 + sp.Icyy*sind(thetad)^2 - sp.Icxy*sind(2*thetad);
Inn = sp.Icxx*sind(thetad)^2 + sp.Icyy*cosd(thetad)^2 + sp.Icxy*sind(2*thetad);
if (Inn-Iee)/maxIcxxIcyy > eps
  % The computed axis was not the major axis, so add 90 degrees of rotation.
  sp.MajorAxisRotDeg = thetad + 90;
else
  sp.MajorAxisRotDeg = thetad;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function runTestCases()

% TEST 1: Formulas from Blevins (2016), Table 1.5.
fprintf('\n')
fprintf('Expected exact values for an annulus:\n')
fprintf('\n')
a    = 50; % Outer radius (mm)
b    = 45; % Inner radius (mm)
xc   = a;
yc   = a;
A    = pi*(a^2-b^2);
Ixx  = 5/4*pi*a^4 - pi*a^2*b^2 - pi/4*b^4;
Iyy  = Ixx;
Icxx = pi/4*(a^4-b^4);
Icyy = Icxx;
Icxy = 0;
deltaThetaDeg = 0.1;
fprintf('Area: %22.15e\n',A)
fprintf('Ixx : %22.15e\n',Ixx)
fprintf('Iyy : %22.15e\n',Iyy)
fprintf('Icxx: %22.15e\n',Icxx)
fprintf('Icyy: %22.15e\n',Icyy)
fprintf('Icxy: %22.15e\n',Icxy)
fprintf('  Xc: %22.15e\n',xc)
fprintf('  Yc: %22.15e\n',yc)
fprintf('\n')
fprintf('Polygon deltaTheta = %f degrees\n',deltaThetaDeg)
fprintf('\n')
thacw = (0:+deltaThetaDeg:+360)';
thcw  = (0:-deltaThetaDeg:-360)';
xyouter = a*[cosd(thacw), sind(thacw)] + a;
xyinner = b*[cosd(thcw ), sind(thcw )] + a;
xy = [xyouter;xyinner;[xyouter(1,1),xyouter(1,2)]];
polygonsp(xy(:,1),xy(:,2))
fprintf('\n')

% TEST 2: Formulas from Blevins (2016), Table 1.5.
fprintf('\n')
fprintf('Expected values for an inclined rectangle:\n')
fprintf('\n')
a      = 50;
b      = 70;
I1     = b*a^3/12;
I2     = a*b^3/12;
thetad = 30;
fprintf('             Area: %22.15e\n',a*b)
fprintf('              Ixx: %22.15e\n',a*b/12*(a^2*sind(thetad)^2+b^2*cosd(thetad)^2))
fprintf('              Iyy: %22.15e\n',a*b/12*(a^2*cosd(thetad)^2+b^2*sind(thetad)^2))
fprintf('              Izz: %22.15e\n',(a^3*b+a*b^3)/12)
fprintf('             Icxy: %22.15e\n',a*b/24*(b^2-a^2)*sind(2*thetad))
fprintf('               Xc: %22.15e\n',0)
fprintf('               Yc: %22.15e\n',0)
fprintf('ImaxPrincipalAxes: %22.15e\n',max(I1,I2))
fprintf('IminPrincipalAxes: %22.15e\n',min(I1,I2))
fprintf('  MajorAxisRotDeg: %22.15e\n',thetad)
fprintf('\n')
xy = [a/2,b/2; -a/2,b/2; -a/2,-b/2; a/2,-b/2; a/2,b/2];
x  = xy(:,1)*cosd(thetad)-xy(:,2)*sind(thetad);
y  = xy(:,1)*sind(thetad)+xy(:,2)*cosd(thetad);
xy = [x(:) y(:)];
polygonsp(xy(:,1),xy(:,2))
fprintf('\n')

% TEST 3: Formulas from Blevins (2016), Table 1.5.
fprintf('\n')
fprintf('Expected values for an equal-length angle section:\n')
fprintf('\n')
thetad = 45;
a   = 4;
b   = 4;
a1  = 1;
b1  = 1;
A   = a1*b+a*b1-a1*b1;
Ixx = 1/3*((a-a1)*b1^3+a1*b^3);
Iyy = 1/3*((b-b1)*a1^3+a^3*b1);
Ixy = 1/4*(a^2*b1^2+a1^2*b^2-a1^2*b1^2);
xc  = (a^2*b1+a1^2*(b-b1))/(2*A);
yc  = (a1*b^2+b1^2*(a-a1))/(2*A);
fprintf('           Area: %22.15e\n',A)
fprintf('            Ixx: %22.15e\n',Ixx)
fprintf('            Iyy: %22.15e\n',Iyy)
fprintf('            Ixy: %22.15e\n',Ixy)
fprintf('           Icxx: %22.15e\n',Ixx-A*yc^2)
fprintf('           Icyy: %22.15e\n',Iyy-A*xc^2)
fprintf('           Icxy: %22.15e\n',Ixy-A*xc*yc)
fprintf('             Xc: %22.15e\n',xc)
fprintf('             Yc: %22.15e\n',yc)
fprintf('MajorAxisRotDeg: %22.15e\n',thetad)
fprintf('\n')
xy = [0,0; a,0; a,b1; a1,b1; a1,b; 0,b; 0,0];
polygonsp(xy(:,1),xy(:,2))
fprintf('\n')

% TEST 4: Formulas from Blevins (2016), Table 1.5.
fprintf('\n')
fprintf('Expected values for an unequal-length angle section:\n')
fprintf('\n')
thetad = 13.8867; % Obtained using an online calculator.
a   = 4;
b   = 8;
a1  = 1;
b1  = 1;
A   = a1*b+a*b1-a1*b1;
Ixx = 1/3*((a-a1)*b1^3+a1*b^3);
Iyy = 1/3*((b-b1)*a1^3+a^3*b1);
Ixy = 1/4*(a^2*b1^2+a1^2*b^2-a1^2*b1^2);
xc  = (a^2*b1+a1^2*(b-b1))/(2*A);
yc  = (a1*b^2+b1^2*(a-a1))/(2*A);
fprintf('           Area: %22.15e\n',A)
fprintf('            Ixx: %22.15e\n',Ixx)
fprintf('            Iyy: %22.15e\n',Iyy)
fprintf('            Ixy: %22.15e\n',Ixy)
fprintf('           Icxx: %22.15e\n',Ixx-A*yc^2)
fprintf('           Icyy: %22.15e\n',Iyy-A*xc^2)
fprintf('           Icxy: %22.15e\n',Ixy-A*xc*yc)
fprintf('             Xc: %22.15e\n',xc)
fprintf('             Yc: %22.15e\n',yc)
fprintf('MajorAxisRotDeg: %22.15e\n',thetad)
fprintf('\n')
xy = [0,0; a,0; a,b1; a1,b1; a1,b; 0,b; 0,0];
polygonsp(xy(:,1),xy(:,2))
fprintf('\n')

end