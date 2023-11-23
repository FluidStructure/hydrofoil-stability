%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate the hydrodynamic integrals for each node

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Nh < 1
    CLnodes = CLv';
else
    CLnodes = [CLv';CLh(2:Nh+1)';CLh(2:Nh+1)'];
end

Lsmat = zeros(Nv+Nh+Nh,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hydrodynamic loading of the vertical section

for i = 1:Nv
    if (i == 1)||(i == Nv)
        Ls = Lv./(2*(Nv-1));
        Lsmat(i,1) = Ls;
    else
        Ls = Lv./(Nv-1);
        Lsmat(i,1) = Ls;
    end
    CL = CLnodes(i);
    % Calculate the centroid of the foil
    [xc,yc] = calcGeom(NACA,HALF_NPANELS,CL);
    % Calculate the hydrodynamic integrals
    [F] = hydroIntegrals(xc,yc,CL,NACA,HALF_NPANELS,rho,Uinf,Ls);
    %F = F.*0;   % Debugging only
    % (For simplicity we will ignore all constant forcing terms)

    % Put these elements into the mass, damping and stiffness matrices
    [Hmass,Hdamp,Hstif] = addHydroIntegralsVertical(Hmass,Hdamp,Hstif,F,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hydrodynamic loading of the horizontal section

if Nh > 0 % First the mid-section (nvth node)
    Ls = (Lh./Nh);   CL = CL1;

    [xc,yc] = calcGeom(NACA,HALF_NPANELS,CL);
    
    [F] = hydroIntegrals(xc,yc,CL,NACA,HALF_NPANELS,rho,Uinf,Ls);
    %F = F.*0;   % Debugging only

    [Hmass,Hdamp,Hstif] = addHydroIntegralsHorizontal(Hmass,Hdamp,Hstif,F,i);
    
    for i = (Nv+1):(Nv+Nh+Nh) % Now for the outlying nodes
        if (i == Nv + Nh)||(i == Nv + Nh + Nh)
            Ls = (Lh./(2*Nh));
            Lsmat(i,1) = Ls;
        else
            Ls = Lh./Nh;
            Lsmat(i,1) = Ls;
        end
        CL = CLnodes(i);
        [xc,yc] = calcGeom(NACA,HALF_NPANELS,CL);
        [F] = hydroIntegrals(xc,yc,CL,NACA,HALF_NPANELS,rho,Uinf,Ls);
        %F = F.*0;   % Debugging only
        [Hmass,Hdamp,Hstif] = addHydroIntegralsHorizontal(Hmass,Hdamp,Hstif,F,i);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
