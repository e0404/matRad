function status = matRad_saveStructs(path, engine)
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

load(path);

if engine=='matlab'
    save('ct.mat', 'ct');
    save('cst.mat', 'cst');
else
    save('ct.mat', '-mat7-binary', 'ct');
    save('cst.mat', '-mat7-binary', 'cst');
end

%Choosing the engine is necessary because Octave has trouble reading .mat files. Doesn't recognize them as binary.

status = 'Files written';

end
