function [readers,writers] = matRad_supportedBinaryFormats()
% matRad function to obtain supported binary formats
% 
% call
%   [read,write] = matRad_supportedBinaryFormats()
%
% input
%
% output
%   read    cell array with file filter in first column, name in second
%           column, and handle to read function in third column
%   write   cell array with file filter in first column, name in second
%           column, and handle to write function in third column
%
% References
%   [1] http://teem.sourceforge.net/nrrd/format.html5
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

%available readers
readers(1).fileFilter = '*.nrrd';
readers(1).name = 'NRRD';
readers(1).handle = @matRad_readNRRD;

readers(2).fileFilter = '*.nii;*.nii.gz';
readers(2).name = 'NifTI';
readers(2).handle = @matRad_readNifTI;

readers(3).fileFilter = '*.mha;*.mhd';
readers(3).name = 'MHA/MHD';
readers(3).handle = @matRad_readMHD;

%available writers
writers(1).fileFilter = '*.nrrd';
writers(1).name = 'NRRD';
writers(1).handle = @matRad_writeNRRD;

writers(2).fileFilter = '*.nii';
writers(2).name = 'NifTI';
writers(2).handle = @matRad_writeNifTI;

writers(3).fileFilter = '*.vtk';
writers(3).name = 'VTK';
writers(3).handle = @matRad_writeVTK;

writers(4).fileFilter = '*.mha';
writers(4).name = 'MHA';
writers(4).handle = @matRad_writeMHA;

writers(5).fileFilter = '*.mhd';
writers(5).name = 'MHD';
writers(5).handle = @matRad_writeMHD;

end

