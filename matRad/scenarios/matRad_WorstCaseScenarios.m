classdef matRad_WorstCaseScenarios < matRad_GriddedScenariosAbstract
%  matRad_WorstCaseScenarios
%  Implements single worst-case shifts per dimension.%  
%
% constructor
%   matRad_WorstCaseScenarios()
%   matRad_WorstCaseScenarios(ct)
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
    

    properties (SetAccess=protected)
        name = 'wcScen';
    end

    properties (Hidden)
        numOfSetupGridPoints = 3;
        numOfRangeGridPoints = 3;
    end
    
    methods
        function this = matRad_WorstCaseScenarios(ct)           
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
        
        function scenarios = updateScenarios(this)
            %Use the static gridded shift function from
            %ImportanceScenarios. We set inclusion of nominal scenarios to
            %false and handle it automatically via the grid point number
            scenarios = this.updateScenarios@matRad_GriddedScenariosAbstract(); 
        end
    end


end
