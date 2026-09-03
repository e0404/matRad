function matRad_makeTallyMCNP(this, ct, fileID_C_rest, binIntervals)
% Description goes here.
%
% call
%   matRad_makeTallyMCNP()
%
% input
%   stf:
%
% output:
%
% References
%   [1] PELOWITZ, D. B., et al. MCNP6 User's Manual. LACP-00634, May, 2013.
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 11/2018
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018-2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matRad_cfg = MatRad_Config.instance();

matRad_cfg.dispInfo('*****\n');
matRad_cfg.dispInfo('Tally type: TMESH3...\n');
matRad_cfg.dispInfo('*****\n');
meshTally.typeCard = 'TMESH\n';
meshTally.geometry = 'RMESH3 %s\n';
meshTally.corA = 'CORA3 %.4f %dI %.4f\n';
meshTally.corB = 'CORB3 %.4f %dI %.4f\n';
meshTally.corC = 'CORC3 %.4f %dI %.4f\n';
meshTally.tallyKeyword = 'TOTAL';

% Write to text file
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, 'C C: Heating tally (one tally located in each voxel of the CT-data)\n');
fprintf(fileID_C_rest, 'C ***************************************************************\n');
fprintf(fileID_C_rest, meshTally.typeCard);
fprintf(fileID_C_rest, meshTally.geometry, meshTally.tallyKeyword);
% Keep in mind matRad LPS coordinate system convention
fprintf(fileID_C_rest, meshTally.corA, .5 * ct.doseGridCT.x_MCNP, ct.doseGridCT.cubeDim(2) - 1, ct.doseGridCT.cubeDim(2) * ct.doseGridCT.x_MCNP + .5 * ct.doseGridCT.x_MCNP);    % Caution: MATLAB indexing
fprintf(fileID_C_rest, meshTally.corB, .5 * ct.doseGridCT.y_MCNP, ct.doseGridCT.cubeDim(1) - 1, ct.doseGridCT.cubeDim(1) * ct.doseGridCT.y_MCNP + .5 * ct.doseGridCT.y_MCNP);
fprintf(fileID_C_rest, meshTally.corC, .5 * ct.doseGridCT.z_MCNP, ct.doseGridCT.cubeDim(3) - 1, ct.doseGridCT.cubeDim(3) * ct.doseGridCT.z_MCNP + .5 * ct.doseGridCT.z_MCNP);

fprintf(fileID_C_rest, 'ENDMD\n');

% Add tallies for RBE weighted dose here

fprintf(fileID_C_rest, 'PRINT 110\n');
fprintf(fileID_C_rest, ['PRDMP ', num2str(ceil(this.config.Num_Primaries)), ' ', num2str(ceil(this.config.Num_Primaries)), ' 1 ', num2str(ceil(this.config.Num_Primaries)), '\n']);  % Control optional MCTAL textfile and set # dumps in RUNTPE to 1

end
