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
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(obj.dicomDir)
    mkdir(obj.dicomDir);
end

matRad_cfg = MatRad_Config.instance();
matRad_cfg.dispInfo('Exporting DICOM for scenario %d to %s:\n',obj.exportScenario,obj.dicomDir);

%Name of Patient & Study
%CT Series
obj = matRad_exportDicomCt(obj);

%RTStruct Series

if ~isempty(obj.cst)
   obj = matRad_exportDicomRTStruct(obj);    
end

%RT Dose Series (Before Plan to have dose reference ids)
if ~isempty(obj.resultGUI)
    obj = matRad_exportDicomRTDoses(obj);
end

%RT Plan at the end
if ~isempty(obj.pln) && ~isempty(obj.resultGUI)
    if obj.enableRtPlanExport
        obj = matRad_exportDicomRTPlan(obj);
    else
        matRad_cfg.dispWarning('pln and resultGUI objects provided, but exporting of RTPlan is disabled by default. Enable it with setting the property enableRtPlanExport to ''true'' before calling the export.');
    end
end

end

