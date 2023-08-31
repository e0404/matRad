function status = matRad_saveStructs(load_path, save_path, engine)
% matRad_saveStructs Mat file transfer for python interface
%
% input
%   path:	storage path
% output
%		d:		status
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2022 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(load_path);

if engine=='matlab'
    save(append(save_path, 'ct.mat'), 'ct');
    save(append(save_path, 'cst.mat'), 'cst');
else
    save(append(save_path, 'ct.mat'), '-mat7-binary', 'ct');
    save(append(save_path, 'ct.mat'), '-mat7-binary', 'cst');
end

%Choosing the engine is necessary because Octave has trouble reading .mat files. Doesn't recognize them as binary.

status = 'Files written';

end
