matRad_rc

ixOAR = 1;
ixPTV = 2;

% define general VOI properties
cst{ixOAR,1} = 0;
cst{ixOAR,2} = 'contour';
cst{ixOAR,3} = 'OAR';

cst{ixPTV,1} = 1;
cst{ixPTV,2} = 'target';
cst{ixPTV,3} = 'TARGET';

% define optimization parameter for both VOIs
cst{ixOAR,5}.TissueClass = 1;
cst{ixOAR,5}.alphaX      = 0.1000;
cst{ixOAR,5}.betaX       = 0.0500;
cst{ixOAR,5}.Priority    = 2;
cst{ixOAR,5}.Visible     = 1;
cst{ixOAR,5}.visibleColor     = [1 0 0];
cst{ixOAR,6}{1,1}.className   = 'DoseObjectives.matRad_SquaredOverdosing';
cst{ixOAR,6}{1,1}.parameters{1}  = 5;
cst{ixOAR,6}{1,1}.penalty     = 100;


cst{ixPTV,5}.TissueClass = 1;
cst{ixPTV,5}.alphaX      = 0.1000;
cst{ixPTV,5}.betaX       = 0.0500;
cst{ixPTV,5}.Priority    = 1;
cst{ixPTV,5}.Visible     = 1;
cst{ixPTV,5}.visibleColor     = [0 1 0];
cst{ixPTV,6}{1,1}.className   = 'DoseObjectives.matRad_SquaredOverdosing';
cst{ixPTV,6}{1,1}.parameters{1}  = 60;
cst{ixPTV,6}{1,1}.penalty     = 800;



 %% Input parameters

%% Create a CT image series
xDim = 700;
yDim = 120;
zDim = 120;

cubeDim      = [xDim yDim zDim];
ct.cube{1} = ones(cubeDim) * 1;
ct.cube{1}(1,1,1) = 0; 

ct.resolution.x = 0.5;
ct.resolution.y = 0.5;
ct.resolution.z = 0.5;

ct.cubeDim = cubeDim;

ct.numOfCtScen  = 1;


%% Lets create either a cubic or a spheric phantom
iso = [600,60,60];

% create an ct image series with zeros - it will be filled later
ct.cubeHU{1} = ones(ct.cubeDim) * 0;
ct.cubeHU{1}(1,1,1) = -1000; 


ct.hlut = [1,0;0,-1024]
% create body of full phantom size
mask = ones(ct.cubeDim);
cst{1,4}{1} = find(mask == 1);

%create target
centerP_corr = iso;
height_corr = 10;
width_corr = 10;
depth_corr = 10;


mask = zeros(ct.cubeDim);

for i=-height_corr/2+1:height_corr/2
    for j=-width_corr/2:width_corr/2
        for k=-depth_corr/2:depth_corr/2
            mask(centerP_corr(1)+i,centerP_corr(2)+j,centerP_corr(3)+k) = 1;
        end
    end
end
cst{2,4}{1} = find(mask == 1);

disp('Done!');
% 
% %%
% figure
% plane = 3;
% slice = iso(3);
% matRad_plotCtSlice(gca,ct.cubeHU,1,plane,slice)
% matRad_plotVoiContourSlice(gca,cst,ct.cube,1,[],plane,slice)
% plot(iso(2),iso(1),'x')

%%
ct1 = ct;
cst1 = cst;

clearvars -except ct cst
save('BOXPHANTOM_LUNG_NARROW_NEW.mat')
