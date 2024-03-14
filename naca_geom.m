clear

if ~exist('getfr.m','file')
    addpath(genpath('lib')); addpath(genpath('data'));
end

NACA = '6512';      % NACA foil shape (4-digit for now)
CL = 1;             % Chord length of the foil (at node points)
HALF_NPANELS = 60;  % Half the number of panels for the NACA foil

iaf.designation=num2str(NACA); iaf.n=HALF_NPANELS;

iaf.HalfCosineSpacing=1; iaf.wantFile=0; iaf.is_finiteTE=0;
iaf.datFilePath='./'; % Current folder

af = naca4gen(iaf);

PANEND(:,1) = flipud(af.x).*CL; PANEND(:,2) = flipud(af.z).*CL;

%fileID = fopen('test.txt','w');
%fprintf(fileID,'%12.8e %12.8e \n',PANEND')
%panend=load('2512.txt');

%pyenv(Version="C:\Users\howellr1\AppData\Local\anaconda3\python")
PANEND(1,2)=0; PANEND(end,2)=0;
[xc,yc,A,Ixx,Iyy,J,x_sc,y_c]=...
    pyrunfile("polygon_sp.py",...
    ["cx","cy","area","ixx_c","iyy_c","j","x_se","y_se"],...
    x=PANEND(:,1),y=PANEND(:,2),materialflg=0);


