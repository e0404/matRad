function engine = getEngineFromPln(pln)
%GETENGINE Summary of this function goes here
%   Detailed explanation goes here

matRad_cfg = MatRad_Config.instance();

engine = [];

initDefaultEngine = false;
%get all available engines for given pln struct, could be done conditional
classList = DoseEngines.matRad_DoseEngineBase.getAvailableEngines(pln);

% Check for a valid engine, and if the given engine isn't valid set boolean 
% to initiliaze default engine at the end of this function
if isfield(pln,'propDoseCalc') && isa(pln.propDoseCalc, 'DoseEngines.matRad_DoseEngineBase')
   engine = pln.propDoseCalc;
elseif isfield(pln,'propDoseCalc') && isstruct(pln.propDoseCalc) && isfield(pln.propDoseCalc,'engine')         
    
    if ischar(pln.propDoseCalc.engine) || isstring(pln.propDoseCalc.engine)
        matchEngines = strcmpi({classList(:).shortName},pln.propDoseCalc.engine);
        if any(matchEngines)
            %instantiate engine
            engineHandle = classList(matchEngines).handle;
            engine = engineHandle(pln);
            
            %engine.assignPropertiesFromPln(pln); %TODO: could this be in the constructor?
        else
            initDefaultEngine = true;
        end
    else
        initDefaultEngine = true;
        matRad_cfg.dispWarning('pln.propDoseCalc.engine field not valid!');
    end
else
    initDefaultEngine = true;
end

% trying to use a default engine which fits
% the given radiation mode, when no valid engine was defined.
% Default Engines are defined in matRad_Config.
if initDefaultEngine
    matchEngines = ismember({classList(:).shortName},matRad_cfg.defaults.propDoseCalc.engine);
    if any(matchEngines)
        engineHandle = classList(matchEngines).handle;

        % unlikely event that multiple engines fit just take the first
        if length(engineHandle) > 1
            engineHandle = engineHandle{1};
        end
        engine = engineHandle(pln);
        matRad_cfg.dispWarning('Using default dose calculation engine %s!', engine.name);
    elseif ~isempty(classList)
        engineHandle = classList(1).handle;
        engine = engineHandle(pln);
        matRad_cfg.dispWarning('Default dose calculation engine not available! Using %s.', engine.name);
    else
        matRad_cfg.dispError('Default dose engine not found!');
    end
end

if isempty(engine)
    matRad_cfg.dispError('No suitable dose engine found!');
end

end

