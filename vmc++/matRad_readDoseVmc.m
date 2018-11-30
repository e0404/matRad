function [bixelDose,bixelDoseError] = matRad_readDoseVmc(filename)
% matRad binary dose import from vmc++
% 
% call
%   [bixelDose,bixelDoseError] = matRad_readDoseVmc(filename)
%
% input
%   filename:   path of input file
%
% output
%   bixelDose       = vector of imported dose values, [D]      = 10^-(10) Gy cm^2
%   bixelDoseError  = vector of imported dose errors, [deltaD] = 10^-(10) Gy cm^2
%
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(filename,'r');

% read header (no regions, no histories, no batches, no beamlets, format specifier (dump_dose))
Header     = fread(fid,5,'int32');
no_regions = Header(1);
dump_dose  = Header(5);

% read dose array
if dump_dose == 2
    dmax            = fread(fid, 1, 'double');
    bixelDose       = fread(fid, no_regions, 'uint16');
    bixelDose       = bixelDose/65534*dmax; % conversion short integers to floating numbers
    bixelDoseError  = 0;
elseif dump_dose == 1
    bixelDose       = fread(fid, no_regions, 'float32');
    bixelDoseError  = fread(fid, no_regions, 'float32');
end

fclose(fid);
return;
