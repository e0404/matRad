function PSI_ijlm =  matRad_calcSecLatMom(vLaSi11,vLaSi22,vLaSi12,vLaSi21,Dev_j,Dev_m) %#codegen
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad function to calculate the second raw moment. Please note that i
% and l depict voxel indices and j and m pencil beam indices. This function
% allows to calculate the the correlation between spot j and multiple
% spots m simultaniously. 
% 
% call
%   PSI_ijlm =  matRad_calcSecLatMom(vLaSi11,vLaSi22,vLaSi12,vLaSi21,Dev_j,Dev_m)
%
% input
%   vLaSi11:        matRads resultGUI struct
%   vLaSi22:        matRad plan meta information struct
%   vLaSi12:        matRad dij dose influence struct
%   vLaSi21:        matRad steering information struct
%   Dev_j:          matRad critical structure struct cst
%   Dev_m

% output
%   PSI_ijlm:         matRad critical structure struct cst
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(isa(vLaSi11,'double'));
assert(isa(vLaSi22,'double'));  
assert(isa(vLaSi12,'double'));  
assert(isa(vLaSi21,'double'));  
assert(isa(Dev_j,'double'));  
assert(isa(Dev_m,'double'));  

coder.varsize('vLaSi11',[1 1]);
coder.varsize('vLaSi22',[Inf 1],[1 0]);
coder.varsize('vLaSi12',[Inf 1],[1 0]);
coder.varsize('vLaSi21',[Inf 1],[1 0]);
coder.varsize('Dev_j',[1 1]);
coder.varsize('Dev_m',[Inf 1],[1 0]);
        
Det           = abs((vLaSi11 .* vLaSi22) - (vLaSi12 .* vLaSi21));
FracDet       = (1./(2*pi.*real(sqrt(Det))));
ExpTerm       = -.5 .* ((Dev_j.*(vLaSi22./Det) + Dev_m.*(-vLaSi12./Det)).* Dev_j +...
                        (Dev_j.*(-vLaSi21./Det) + Dev_m.*(vLaSi11./Det)).* Dev_m);
PSI_ijlm      = FracDet .* exp(ExpTerm);
                                
end

