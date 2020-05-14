function topasDose = matRad_getWslTopasData(wslDistribution,userName)

load(['\\wsl$\',wslDistribution,'\home\',userName,'\matfiles\MCexport\MCparam.mat']);
topasDose.MCparam = MCparam;

topasSum = zeros(MCparam.cubeDim(2),MCparam.cubeDim(1),MCparam.cubeDim(3));
for k = 1:MCparam.nbRuns
    data = csvread(['\\wsl$\',wslDistribution,'\home\',userName,'\matfiles\MCexport\simData_matRad_plan_field1_run',num2str(k),'_physicalDose.csv'],8,0);
    for n = 1:length(data)
        dose{k}(data(n,1)+1,data(n,2)+1,data(n,3)+1) = data(n,4);
    end
    topasSum = topasSum + dose{k};
end
correctionFactor = double(MCparam.nbParticles) / double(MCparam.nbHistories);

% correct the topas calculated dose
topasDose.physicalDose = zeros(MCparam.cubeDim(2),MCparam.cubeDim(1),MCparam.cubeDim(3));
topasDose.physicalDose(:) = correctionFactor(:).*topasSum(:);

% switch coordinate systems from topas (x,y,z) to matlab (y,x,z)
topasDose.physicalDose = permute(topasDose.physicalDose,[2 1 3]);

end