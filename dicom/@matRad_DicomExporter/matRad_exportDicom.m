function obj = matRad_exportDicom(obj)
% matRad function to export current workspace to DICOM. 
% Function of matRad_DicomExporter
% 
% call
%   matRad_DicomExporter.matRad_exportDicom()
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(obj.dicomDir)
    mkdir(obj.dicomDir);
end

%Name of Patient & Study
%CT Series
obj = matRad_exportDicomCt(obj);

%RTStruct Series

if ~isempty(obj.cst)
   obj = matRad_exportDicomRTStruct(obj);    
end

%if ~isempty(obj.pln)
%    obj = matRad_exportDicomRTPlan(obj);
%end

if ~isempty(obj.resultGUI)
    obj = matRad_exportDicomRTDoses(obj);
end

%Pln Series

%Dose Series
end

