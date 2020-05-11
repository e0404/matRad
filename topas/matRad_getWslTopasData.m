function topasDose = matRad_getWslTopasData(wslDistribution,userName,ct)

load(['\\wsl$\',wslDistribution,'\home\',userName,'\matfiles\MCexport\MCparam.mat']);
topasDose.MCparam = MCparam;

topasSum = zeros(MCparam.cubeDim(2),MCparam.cubeDim(1),MCparam.cubeDim(3));
for k = 1:MCparam.nbRuns
    data = csvread(['\\wsl$\',wslDistribution,'\home\',userName,'\matfiles\MCexport\simData_matRad_plan_field1_run',num2str(k),'_physicalDose.csv'],8,0);
    for n = 1:length(data)
        dose{k}(data(n,1)+1,data(n,2)+1,data(n,3)+1) = data(n,4);
    end
    topasSum = topasSum + dose{k}/MCparam.nbRuns;
end
voxelVolume = 1.0e-3*(MCparam.voxelDimensions.x*MCparam.voxelDimensions.y*MCparam.voxelDimensions.z);
correctionFactor = double(MCparam.nbParticles) / double(MCparam.nbHistories);

topasDose.physicalDose = zeros(MCparam.cubeDim(2),MCparam.cubeDim(1),MCparam.cubeDim(3));
topasDose.physicalDose(:) = correctionFactor(:).*topasSum(:);

topasDose.physicalDose = rot90(topasDose.physicalDose,-1);
topasDose.physicalDose = flip(topasDose.physicalDose,2);

end