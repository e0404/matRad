function [bixelDose,bixelDoseError] = matRad_readDoseVmc(filename,VmcOptions)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(filename,'r');

% read header (no regions, no histories, no batches, no beamlets, format specifier (dump_dose))
switch VmcOptions.run.version
    case 'Carleton'
        Header      = fread(fid,1,'int32');
        no_regions  = Header(1);
        dump_dose   = VmcOptions.scoringOptions.outputOptions.dumpDose;
    case 'dkfz'
        Header      = fread(fid,5,'int32');
        no_regions  = Header(1);
        dump_dose   = Header(5);
end

% read dose array
if dump_dose == 2
    dmax            = fread(fid, 1, 'double');
    bixelDose       = fread(fid, no_regions, 'uint16');
    bixelDose       = bixelDose/65534*dmax; % conversion short integers to floating numbers
    bixelDoseError  = zeros(size(bixelDose));
elseif dump_dose == 1 
    bixelDose       = fread(fid, no_regions, 'float32');
    bixelDoseError  = fread(fid, no_regions, 'float32');
end
fclose(fid);

% reshape into array, permute y <-> x, reshape back into column
bixelDose = reshape(bixelDose,VmcOptions.geometry.dimensions([2 1 3]));
bixelDose = permute(bixelDose,[2 1 3]);
bixelDose = reshape(bixelDose,[],1);

bixelDoseError = reshape(bixelDoseError,VmcOptions.geometry.dimensions([2 1 3]));
bixelDoseError = permute(bixelDoseError,[2 1 3]);
bixelDoseError = reshape(bixelDoseError,[],1);

end




