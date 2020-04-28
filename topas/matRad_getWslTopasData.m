function topasDose = matRad_getWslTopasData(wslDistribution,userName,TopasConfig)

param = load(['\\wsl$\',wslDistribution,'\home\',userName,'\matfiles\MCexport\MCparam.mat']);
topasDose.MCparam = param.MCparam;

for k = 1:param.MCparam.nbRuns
data = csvread([TopasConfig.filepath,'simData_matRad_plan_field1_run',num2str(k),'_physicalDose.csv'],8,0);
dose{k} = data(:,4);
end

topasDose = sum(cell2mat(dose),2);
topasDose = reshape(topasDose(:),160,100,100);

end