function classList = getAvailableEngines(pln,optionalPaths)
% Returns a list of names and coresponding handle for available dose calc engines
%   Returns all dose calc engines in the package when no arg is
%   given. If no engines are found return gonna be empty.
%
% call:
%   classList = DoseEngines.matRad_DoseEngine.getAvailableEngines(pln,optional_path)
%
% input:
%   pln:            containing proposed dose calc and machine file informations
%   optionalPath:   cell array of other folders to search in
%
% returns:
%   classList: struct array containing name, shortName, className, and
%   handle for construction of the Engine

matRad_cfg = MatRad_Config.instance();

%Parse inputs
if nargin < 2
    optionalPaths = {};
else
    if ~(iscellstr(optionalPaths) && all(optionalPaths))
        matRad_cfg.dispError('Invalid path array!');
    end
end

if nargin < 1
    pln = [];
else
    if ~(isstruct(pln) || isempty(pln))
        matRad_cfg.dispError('Invalid pln!');
    end
end


%Get available, valid classes through call to matRad helper function
%for finding subclasses
availableDoseEngines = matRad_findSubclasses('DoseEngines.matRad_DoseEngineBase','packages',{'DoseEngines'},'folders',optionalPaths,'includeAbstract',false);

%Now filter for pln
ix = [];

if nargin >= 1 && ~isempty(pln)
    machine = matRad_loadMachine(pln);
    machineMode = machine.meta.radiationMode;

    for cIx = 1:length(availableDoseEngines)
        mc = availableDoseEngines{cIx};
        availabilityFuncStr = [mc.Name '.isAvailable'];
        %availabilityFunc = str2func(availabilityFuncStr); %str2func  does not seem to work on static class functions in Octave 5.2.0
        try
            %available = availabilityFunc(pln,machine);
            available = eval([availabilityFuncStr '(pln,machine)']);
        catch
            available = false;
            mpList = mc.PropertyList;
            if matRad_cfg.isMatlab
                loc = find(arrayfun(@(x) strcmp('possibleRadiationModes',x.Name),mpList));
                propValue = mpList(loc).DefaultValue;
            else
                loc = find(cellfun(@(x) strcmp('possibleRadiationModes',x.Name),mpList));
                propValue = mpList{loc}.DefaultValue;
            end
    
            if any(strcmp(propValue, pln.radiationMode))
                % get radiation mode from the in pln proposed basedata machine file
                % add current class to return lists if the
                % radiation mode is compatible
                if(any(strcmp(propValue, machineMode)))
                    available = true;
                    
                end
            end
        end
        if available
            ix = [ix cIx];
        end
    end

    availableDoseEngines = availableDoseEngines(ix);
end

classNameList = cellfun(@(mc) mc.Name,availableDoseEngines,'UniformOutput',false);

if matRad_cfg.isMatlab
    shortNameList = cellfun(@(mc) mc.PropertyList(strcmp({mc.PropertyList.Name}, 'shortName')).DefaultValue,availableDoseEngines,'UniformOutput',false);
    nameList = cellfun(@(mc) mc.PropertyList(strcmp({mc.PropertyList.Name}, 'name')).DefaultValue,availableDoseEngines,'UniformOutput',false);
else
    %Indexing a cell array with . not possible (in Octave)
    shortNameList = cellfun(@(mc) mc.PropertyList{find(cellfun(@(p) strcmp(p.Name, 'shortName'),mc.PropertyList))}.DefaultValue,availableDoseEngines,'UniformOutput',false);
    nameList = cellfun(@(mc) mc.PropertyList{find(cellfun(@(p) strcmp(p.Name, 'name'),mc.PropertyList))}.DefaultValue,availableDoseEngines,'UniformOutput',false);
end

%make sure the default engines are the first ones listed
for defaultEngine = matRad_cfg.propDoseCalc.defaultDoseEngines
    findDefaultIx = strcmp(defaultEngine,shortNameList);
    if ~isempty(findDefaultIx)
        shortNameList = [shortNameList(findDefaultIx), shortNameList(~findDefaultIx)];
        classNameList = [classNameList(findDefaultIx), classNameList(~findDefaultIx)];
    end
end


handleList = cellfun(@(namestr) str2func(namestr),classNameList,'UniformOutput',false);

classList = cell2struct([shortNameList; nameList; classNameList; handleList],{'shortName','name','className','handle'});

end


