function [xc,yc,Ixx,Iyy,J,A] = calcGeom(NACA,HALF_NPANELS,CL)

% Get NACA section using function naca4gen
iaf.designation=num2str(NACA); iaf.n=HALF_NPANELS; iaf.HalfCosineSpacing=1;
iaf.wantFile=0; iaf.datFilePath='./'; % Current folder
iaf.is_finiteTE=0;

af = naca4gen(iaf);

PANEND(:,1) = flipud(af.x).*CL;
PANEND(:,2) = flipud(af.z).*CL;
PANEND(1,2)=0; PANEND(end,2)=0;
%fileID = fopen('test.txt','w');
%fprintf(fileID,'%12.8e %12.8e \n',PANEND')

% Get the geometry properties
methodflg=1;
if methodflg==1
    [ geom, iner, cpmo ] = polygeom(PANEND(:,1),PANEND(:,2));

    xc = geom(2); yc = geom(3); A = geom(1);
    Ixx = iner(4); Iyy = iner(5); J = cpmo(5);
else
    [xc,yc,A,Ixx,Iyy,J,x_sc,y_c]=...
    pyrunfile("polygon_sp.py",...
    ["cx","cy","area","ixx_c","iyy_c","j","x_se","y_se"],...
    x=PANEND(:,1),y=PANEND(:,2),materialflg=0);
end

keyboard

