%% create water phantom with slab
matRad_rc

% first vertical, second horizontal, third longitudinal
resolution          = [1, 1, 1];
cubeDim             = [150, 50, 50]; 
slabThickness       = 30;  %[mm] 
relStoppingpower    = 1.8; %[mm]
distToEntrance      = 5;  %80 %[mm]


ct.cubeDim = cubeDim;
ct.resolution.x = resolution(1);
ct.resolution.y = resolution(2);
ct.resolution.z = resolution(3);
ct.numOfCtScen = 1;
ct.cube{1} = ones(cubeDim);

horStart = distToEntrance/resolution(1);
horEnd   = (distToEntrance + slabThickness)/resolution(1);

%% create slab reaching half into ray
% for i = horStart:horEnd;
%     for j = 1:round(cubeDim(2)/2)
%         for k = 1:cubeDim(3)
%             ct.cube{1}(i,j,k) = relStoppingpower;
%         end
%     end
% end

%% create diagonal wedge in beam
for i = horStart:horEnd;
    for j = 1:round(cubeDim(2)*(i - horStart) / (horEnd - horStart))
        for k = 1:cubeDim(3)
            ct.cube{1}(i,j,k) = relStoppingpower;
        end
    end
end


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


mask = ones(ct.cubeDim);
cst{1,4}{1} = find(mask == 1);
for i = horStart:horEnd
    for j = round(cubeDim(2)*0.25):round(cubeDim(2)*0.75)
        for k = round(cubeDim(3)*0.25):round(cubeDim(3)*0.75)
            mask(i,j,k) = 2;
        end
    end
end
cst{2,4}{1} = find(mask == 2);

clearvars -except ct cst
matRadGUI