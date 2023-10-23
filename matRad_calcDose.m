function dij = matRad_calcDose(ct,cst,stf,pln)
% matRad dose calculation automaticly creating the appropriate dose engine
% for the given pln struct and called the associated dose calculation funtion
%
% call
%   dij =  matRad_calcDose(ct,stf,pln,cst)
%
% input
%   ct:         ct cube
%   cst:        matRad cst struct
%   stf:        matRad steering information struct
%   pln:        matRad plan meta information struct
%
%
% output
%   dij:            matRad dij struct
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
initDefaultEngine = false;

%Deprecation warnings
if ~isfield(stf,'machine')
    matRad_cfg.dispDeprecationWarning('stf should contain the machine name in the ''machine'' field since matRad 3. Manually adding ''%s'' from pln.',pln.machine);
    for i=1:numel(stf)
        stf(i).machine = pln.machine;
    end
end

%get all available engines for given pln struct, could be done conditional
[nameList, ~, handleList] = DoseEngines.matRad_DoseEngineBase.getAvailableEngines(pln);

if(isfield(pln,'propDoseCalc') && isfield(pln.propDoseCalc,'engine'))


    % check if there is an engine field in the pln struct
    % check if it is a dose calc engine and not some other oject
    if isa(pln.propDoseCalc.engine, 'DoseEngines.matRad_DoseEngineBase')
        % set engine
        engine = pln.propDoseCalc.engine;
    elseif  ischar(pln.propDoseCalc.engine)
        if any(strcmpi(nameList,pln.propDoseCalc.engine))
            %instantiate engine
            engineHandle = handleList{strcmpi(nameList,pln.propDoseCalc.engine)};
            engine = engineHandle(pln);
            
            %engine.assignPropertiesFromPln(pln); %TODO: could this be in the constructor?
        else
            % if the given engine isn't valid set boolean to initiliaze
            % it at the end of this function
            initDefaultEngine = true;
        end
    else
        initDefaultEngine = true;
        matRad_cfg.dispWarning('pln.propDoseCalc.engine field not valid!');
    end
else
    % if the given engine isn't valid set boolean to initiliaze
    % it at the end of this function
    initDefaultEngine = true;
end

% trying to use a default engine which fits
% the given radiation mode, when no valid engine was defined.
% Default Engines are defined in matRad_Config.
if initDefaultEngine

    if any(ismember(nameList,matRad_cfg.propDoseCalc.defaultDoseEngines))
        engineHandle = handleList{ismember(nameList,matRad_cfg.propDoseCalc.defaultDoseEngines)};

        % unlikely event that multiple engines fit just take the first
        if length(engineHandle) > 1
            engineHandle = engineHandle{1};
        end
        engine = engineHandle(pln);
        matRad_cfg.dispWarning('Using default dose calculation engine %s!', engine.name);
    elseif ~isempty(nameList)
        engineHandle = handleList{1};
        engine = engineHandle(pln);
        matRad_cfg.dispWarning('Default dose calculation engine not available! Using %s.', engine.name);
    else
        matRad_cfg.dispError('No dose engine found!');
    end


end
%call the calcDose funktion
dij = engine.calcDose(ct,cst,stf,pln);

end
