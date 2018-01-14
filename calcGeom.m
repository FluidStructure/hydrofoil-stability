function [xc,yc,Ixx,Iyy,J,A] = calcGeom(NACA,HALF_NPANELS,CL)

% Get NACA section using function naca4gen
iaf.designation=num2str(NACA);
iaf.n=HALF_NPANELS;
iaf.HalfCosineSpacing=1;
iaf.wantFile=0;
iaf.datFilePath='./'; % Current folder
iaf.is_finiteTE=0;
af = naca4gen(iaf);
PANEND(:,1) = flipud(af.x).*CL;
PANEND(:,2) = flipud(af.z).*CL;

% Get the geometry properties
[ geom, iner, cpmo ] = polygeom(PANEND(:,1),PANEND(:,2));
xc = geom(2);
yc = geom(3);

A = geom(1);

% Get the moments of area about centroidal x and y axes
Ixx = iner(4);
Iyy = iner(5);
J = cpmo(5);