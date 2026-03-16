%% Calculate opening angles for maximum field size 
% MLC 1

lengthBeamline1 = [62 274 288 428 478 493];
hight1 = [16 18 30 20 20 22];
width1 = [18 24 35 30 30 30];

sourceWidht = 15;   % square dimension of converter plate

funA = @(ho) ho/2 - sourceWidht/2;
funB = @(br) br/2 + sourceWidht/2;
funC = @(ho) ho/2 + sourceWidht/2;
funD = @(br) br/2 - sourceWidht/2;

funU = @(a,b) sqrt(a^2 + b^2);
funV = @(c,b) sqrt(c^2 + b^2);
funW = @(c,d) sqrt(c^2 + d^2);
funX = @(a,d) sqrt(a^2 + d^2);

funTheta = @(ho, br,l) [atand(funX(funA(ho), funD(br))/l), ...
    atand(funU(funA(ho), funB(br))/l), ...
    atand(funV(funB(br), funC(ho))/l), ...
    atand(funW(funA(ho), funD(br))/l)];
    

maxOpeningAngle_MLC1 = 0;
openingAngles_MLC1 = zeros(4,length(lengthBeamline1));
for counter = 1:length(lengthBeamline1)
    dummy = max(max(funTheta(hight1(counter), width1(counter), lengthBeamline1(counter))));
    maxOpeningAngle_MLC1 = max(maxOpeningAngle_MLC1, dummy);
    openingAngles_MLC1(:,counter) = funTheta(hight1(counter), width1(counter), lengthBeamline1(counter))';
end

MLC_1.maxOpeningAngle_MLC1 = maxOpeningAngle_MLC1;
MLC_1.openingAngles_MLC1 = openingAngles_MLC1;

MLC_1.lengthBeamline1 = lengthBeamline1;
MLC_1.hight1 = hight1;
MLC_1.width1 = width1;

%% Calculate opening angles for maximum field size 
% MLC 3

lengthBeamline3 = [62 274 288 463 493 523];
hight3 = [16 18 30 24 18 6]; %18];
width3 = [18 24 35 30 26.5 6]; %26.5];

sourceWidht = 15;   % square dimension of converter plate

funA = @(ho) ho/2 - sourceWidht/2;
funB = @(br) br/2 + sourceWidht/2;
funC = @(ho) ho/2 + sourceWidht/2;
funD = @(br) br/2 - sourceWidht/2;

funU = @(a,b) sqrt(a^2 + b^2);
funV = @(c,b) sqrt(c^2 + b^2);
funW = @(c,d) sqrt(c^2 + d^2);
funX = @(a,d) sqrt(a^2 + d^2);

funTheta = @(ho, br,l) [atand(funX(funA(ho), funD(br))/l), ...
    atand(funU(funA(ho), funB(br))/l), ...
    atand(funV(funB(br), funC(ho))/l), ...
    atand(funW(funA(ho), funD(br))/l)];
    

maxOpeningAngle_MLC3 = 0;
openingAngles_MLC3 = zeros(4,length(lengthBeamline3));
for counter = 1:length(lengthBeamline3)
    dummy = max(max(funTheta(hight3(counter), width3(counter), lengthBeamline3(counter))));
    maxOpeningAngle_MLC3 = max(maxOpeningAngle_MLC3, dummy);
    openingAngles_MLC3(:,counter) = funTheta(hight3(counter), width3(counter), lengthBeamline3(counter))';
end

MLC_3.maxOpeningAngle_MLC3 = maxOpeningAngle_MLC3;
MLC_3.openingAngles_MLC3 = openingAngles_MLC3;

MLC_3.lengthBeamline3 = lengthBeamline3;
MLC_3.hight3 = hight3;
MLC_3.width3 = width3;

%% Clean up 
clearvars -except MLC_1 MLC_3