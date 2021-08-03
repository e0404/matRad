%% conduct gamma analysis
addpath(fileparts(pwd));
matRad_rc

% load basedata 
    load brachy_HDR
    HDRmachine = machine;
    load brachy_LDR
    LDRmachine = machine;
    clear machine

%% 3D coordinates
    % set resolution
    res = 3;
    maxRadius = 120;
    DoseCutoff = 15;
    calcCode = 'here';
    %calcCode = 'ref';
    dim = '2D';
    

    % prepare grid
    x = [-maxRadius:res:maxRadius];
    y = [-maxRadius:res:maxRadius];
    z = 0;

    [grids.x,grids.y,grids.z] = meshgrid(x,y,z);
    grids.dist = sqrt(grids.x.^2 + grids.y.^2 + grids.z.^2);
    theta = matRad_getThetaMatrix([1,0,0],grids);
clear res maxRadius x y z
%% calculate doses
    switch calcCode
        case 'here'
            if strcmp(dim,'2D')
                DoseRate = matRad_getDoseRate2D_poly(HDRmachine,grids.dist,theta);
            else
                DoseRate = matRad_getDoseRate1D_poly(LDRmachine,grids.dist);
            end
        case 'ref'
            if strcmp(dim,'2D')
                HDRmachine.data.AnisotropyMatrix = ...
                HDRmachine.data.AnisotropyFunctionValue; % quick fix
                source = Source(HDRmachine.data);
                DoseRate = source.getDoseRate2D(grids.dist/10,theta);
            else
                %LDRmachine.data.AnisotropyFunctionValue = ...
                %LDRmachine.data.AnisotropyFunctionValue; % quick fix
                source = Source(LDRmachine.data);
                DoseRate = source.getDoseRate1D(grids.dist/10);  
            end
    end
    
    DoseRate(DoseRate>DoseCutoff) = DoseCutoff;
    DoseRate(DoseRate<0) = 0;
    
    refDos.TG43_2D.fullDose = DoseRate;

%% visualize stuff if needed
%     surf(grids.x,grids.y,DoseRate,'LineStyle','none')
%     xlabel('x[mm]')
%     ylabel('y[mm]')
%     zlabel('2D approx DoseRate')
    image(DoseRate,'CDataMapping','scaled')
    colorbar