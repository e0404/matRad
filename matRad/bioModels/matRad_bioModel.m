function model = matRad_bioModel(radiationMode, model, providedQuantities)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  matRad_bioModel
%  This is a helper function to instantiate a matRad_BiologicalModel. This
%  function currently exists for downwards compatability, as the new
%  Biological Models will follow a polymorphic software architecture
%
% call
%   matRad_bioModel(radiationMode, model)
%
%   e.g. pln.bioModel = matRad_bioModel('protons','MCN')
%
% input
%   radiationMode:      radiation modality 'photons' 'protons' 'helium' 'carbon' 'brachy'
%   
%   model:              string to denote which biological model is used
%                       'none': for photons, protons, carbon                'constRBE': constant RBE for photons and protons
%                       'MCN': McNamara-variable RBE model for protons      'WED': Wedenberg-variable RBE model for protons
%                       'LEM': Local Effect Model for carbon ions
%
%   providedQuantities: optional cell string of provided quantities to
%                       check if the model can be evaluated
%
% output
%   model:              instance of a biological model
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    model = matRad_BiologicalModel.validate(model,radiationMode);
else
    model = matRad_BiologicalModel.validate(model,radiationMode,providedQuantities);
end

end % end function