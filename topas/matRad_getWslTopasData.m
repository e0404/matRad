function topasDose = matRad_getWslTopasData(wslDistribution,userName)

topasDose = load(['\\wsl$\',wslDistribution,'\home\',userName,'\matfiles\MCexport\MCdata.mat']);
topasDose = topasDose.MCdata;
param = load(['\\wsl$\',wslDistribution,'\home\',userName,'\matfiles\MCexport\MCparam.mat']);
topasDose.MCparam = param;

end