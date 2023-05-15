function model = matRad_bioModel(sRadiationMode,sQuantityOpt, sModel)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  matRad_bioModel
%  This is a helper function to instantiate a matRad_BiologicalModel. This
%  function currently exists for downwards compatability, as the new
%  Biological Models will follow a polymorphic software architecture
%
% call
%   matRad_bioModel(sRadiationMode,sQuantityOpt, sModel)
%
%   e.g. pln.bioParam = matRad_bioModel('protons','constRBE','RBExD')
%
% input
%   sRadiationMode:     radiation modality 'photons' 'protons' 'carbon'
%   sQuantityOpt:       string to denote the quantity used for
%                       optimization 'physicalDose', 'RBExD', 'effect'
%   sModel:             string to denote which biological model is used
%                       'none': for photons, protons, carbon                'constRBE': constant RBE for photons and protons
%                       'MCN': McNamara-variable RBE model for protons      'WED': Wedenberg-variable RBE model for protons
%                       'LEM': Local Effect Model for carbon ions
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

model = matRad_BiologicalModel(sRadiationMode,sQuantityOpt,sModel);

end % end class definition