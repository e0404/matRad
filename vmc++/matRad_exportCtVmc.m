function matRad_exportCtVmc(ct,filename)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad binary CT export for vmc++
% 
% call
%   matRad_exportCtVmc(ct,filename)
%
% input
%   ct:             matRad ct struct
%   filename:       path where CTfile is created
%
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(filename,'w');

% write ct dimensions
fwrite(fid,ct.cubeDim([2 1 3]),'int32');

% write voxel corner location in cm in physical cs with ct cube corner at [.5 .5 .5]
X = (0.5:(ct.cubeDim(2)+0.5))*ct.resolution.x/10;
Y = (0.5:(ct.cubeDim(1)+0.5))*ct.resolution.y/10;
Z = (0.5:(ct.cubeDim(3)+0.5))*ct.resolution.z/10;

fwrite(fid,X,'float32');
fwrite(fid,Y,'float32');
fwrite(fid,Z,'float32');

% write voxel densities
% first permute indices y <-> x
ctVMC = permute(ct.cube{1},[2 1 3]);
% then reshape into single column
ctVMC = reshape(ctVMC,[],1);
fwrite(fid,ctVMC,'float32');

fclose(fid);
