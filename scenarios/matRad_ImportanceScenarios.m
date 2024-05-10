classdef matRad_ImportanceScenarios < matRad_GriddedScenariosAbstract
%  matRad_ImportanceScenarios
%  Implements gridded importance scenarios, i.e., weighted according to a
%  probability distribution. It is not advised to create "combined" grids
%  with a large number of grid-points as the curse of dimensionality will
%  quickly break memory requirements when putting this in to dose influence
%  matrix computation.
%  
%
% constructor
%   matRad_ImportanceScenarios()
%   matRad_ImportanceScenarios(ct)
%
% input
%   ct:                 ct cube
%
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
    
    properties (AbortSet = true)
        numOfSetupGridPoints = 9;
        numOfRangeGridPoints = 9;
    end

    properties (SetAccess=protected)
        name = 'impScen';
    end
    
    methods
        function this = matRad_ImportanceScenarios(ct)           
            if nargin == 0 
                superclassArgs = {};
            else
                superclassArgs = {ct};
            end
            
            this@matRad_GriddedScenariosAbstract(superclassArgs{:});

            %TODO: We could do this automatically in the superclass
            %Octave 5 has a bug there and throws an error
            this.updateScenarios();
        end

        function set.numOfSetupGridPoints(this,numGridPoints)
            valid = isscalar(numGridPoints) && numGridPoints > 0;
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid number of setup grid points, needs to be a positive scalar!');
            end
            this.numOfSetupGridPoints = inumGridPoints;
            this.updateScenarios();
        end

        function set.numOfRangeGridPoints(this,numGridPoints)
            valid = isscalar(numGridPoints) && numGridPoints > 0;
            if ~valid 
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid number of range grid points, needs to be a positive scalar!');
            end
            this.numOfRAngeGridPoints = inumGridPoints;
            this.updateScenarios();
        end
    end
end

