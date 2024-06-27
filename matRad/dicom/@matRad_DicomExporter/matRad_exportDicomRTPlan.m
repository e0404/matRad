function obj = matRad_exportDicomRTPlan(obj)
% matRad function to export resultGUI to dicom.
% 
% call
%   matRad_exportDicomRTDoses(resultGUI,ct,pln,fieldnames)
%
% input
%   resultGUI:      matRad resultGUI struct with different beams. Note that
%                   the summation (called plan) of the beams is named 
%                   without subscripts, e.g. physical_Dose.
%   ct:             matRad ct struct
%
% output
%   resultGUI:      matRad resultGUI struct with different beams. Note that
%                   the summation (called plan) of the beams is named 
%                   without subscripts, e.g. physical_Dose.
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

matRad_cfg = MatRad_Config.instance();
matRad_cfg.dispWarning('RTPlan export is not yet implemented...\n');            

end
