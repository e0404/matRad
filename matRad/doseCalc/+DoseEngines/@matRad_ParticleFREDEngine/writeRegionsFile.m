function writeRegionsFile(this,fName)
% FRED helper to write file for geometrical components
% call
%   writeRegionsFile(fName)
% 
% input
%   fName: string specifying the file path and name for saving the data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matRad_cfg = MatRad_Config.instance();

fID = fopen(fName, 'w');
try
    % Write patient component
    fprintf(fID,'region<\n');
    fprintf(fID,'\tID=Phantom\n');
    fprintf(fID,'\tCTscan=inp/regions/%s\n', this.patientFilename);
    fprintf(fID,'\tO=[%i,%i,%i]\n', 0,0,0);
    fprintf(fID,'\tpivot=[0.5,0.5,0.5]\n');

    % l=e1; u=e2;
    % x in Room coordinates is  x in patient frame
    % y in Romm coordinates is -y in patient frame
    % Voxels in y-direection in matRad grow in -y direction in FRED Room reference 
    fprintf(fID, '\tl=[%1.1f,%1.1f,%1.1f]\n', 1,0,0);
    fprintf(fID, '\tu=[%1.1f,%1.1f,%1.1f]\n', 0,-1,0);

    % Syntax changes for scorers according to direct or ij calculation
    if this.calcDoseDirect
        fprintf(fID,'\tscore=[');
    else
        fprintf(fID,'\tscoreij=[');
    end

    if numel(this.scorers)>1
        for k=1:size(this.scorers,2)-1
            fprintf(fID,'%s,', this.scorers{k});
        end
    end
    fprintf(fID,'%s]\n', this.scorers{end});    
    
    fprintf(fID,'region>\n');

    % Write Room parameters
    fprintf(fID, 'region<\n');
    fprintf(fID, '\tID=Room\n');
    fprintf(fID, '\tmaterial=%s\n', this.roomMaterial);
    fprintf(fID, 'region>\n');
    
    % Write HU table if needed
    switch this.HUtable
        case 'internal'
            fprintf(fID, 'lUseInternalHU2Mat=t\n');
        otherwise
            fprintf(fID, 'include: inp/regions/hLut.inp\n');
            this.writeHlut(this.HUtable);
    end
    
    % Toogle HU clamping if requested
    if this.HUclamping
        fprintf(fID, 'lAllowHUClamping=t\n');
    end

    if ~isempty(this.dijFormatVersion) && this.isVersionHigher('3.70.0')
        fprintf(fID, 'ijFormatVersion = %s\n', this.dijFormatVersion);
    end

catch ME
    matRad_cfg.dispError(['Failed to write regions file. Exit with error: ', ME.message]);

end

fclose(fID);
end