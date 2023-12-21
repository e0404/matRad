function LETd = matRad_calcLETd(dij,w)
% Calculates LETd 
%
% call
%   dij = matRad_calcLETd(dij,w)
%
% input
%   dij:  matRad dij struct
%   w:    weight vector
%
% output
%   dij:  matRad dij struct with LETd
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LETd = zeros(size(dij.physicalDose{1}));
d = dij.physicalDose{1} * w;
wLET = dij.mLETDose{1} * w;
LETd(d>0) = wLET(d>0)./d(d>0);

end