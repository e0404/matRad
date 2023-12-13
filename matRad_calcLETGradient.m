function LETgrad = matRad_calcLETGradient(dij,doseGrad,scen)
% Calculates LET gradient
%
% call
%   LETgrad = matRad_calcLETGradient(dij,doseGrad,scen)
%
% input
%   dij:        matRad dij struct
%   doseGrad:   dose gradient
%   scen:       scen for projection
%
% output
%   LETgrad:   LET gradient
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LETgrad = (((doseGrad{scen}'*dij.physicalDose)' * dij.mLET) * dij.physicalDose ...
            - (dij.mLETDose * (doseGrad{scen}'*dij.physicalDose)')) / ((dij.physicalDose)^2);

end