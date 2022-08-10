function [nameList, classList, handleList] = getAvailableEngines(pln,optionalPaths)
% Returns a list of names and coresponding handle for available dose calc engines
%   Returns all dose calc engines in the package when no arg is
%   given. If no engines are found return gonna be empty.
%
% call:
%   [nameList, handleList] = DoseEngines.matRad_DoseEngine.getAvailableEngines(pln,optional_path)
%
% input:
%   pln:            containing proposed dose calc and machine file informations
%   optionalPath:   cell array of other folders to search in
%
% returns:
%   nameList: cell-array conatining readable names for engines
%   classList: cell-array conatining full classnamens for available engines
%   handleList: cell-array containing function-handles to
%                  available engines constructor (call the included handle by adding Parentheses e.g. handleList{1}())

matRad_cfg = MatRad_Config.instance();

%Parse inputs
p = inputParser;
p.addOptional('pln',[],@(x) isstruct(x) || isempty(x));
p.addOptional('optionalPaths',{},@(x) iscellstr(x) && all(isfolder(x)));
if nargin < 1
    p.parse();
elseif nargin < 2
    p.parse(pln);
else
    p.parse(pln,optionalPaths);
end
pln = p.Results.pln;
optionalPaths = p.Results.optionalPaths;

%Get available, valid classes through call to matRad helper function
%for finding subclasses
availableDoseEngines = matRad_findSubclasses('DoseEngines.matRad_DoseEngine','packages',{'DoseEngines'},'folders',optionalPaths,'includeAbstract',false);

%Now filter for pln
ix = [];

if nargin >= 1 && ~isempty(pln)
    machine = matRad_loadMachine(pln);
    machineMode = machine.meta.radiationMode;

    for cIx = 1:length(availableDoseEngines)
        mc = availableDoseEngines{cIx};
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
                ix = [ix cIx];
            end
        end
    end

    availableDoseEngines = availableDoseEngines(ix);
end

classList = cellfun(@(mc) mc.Name,availableDoseEngines,'UniformOutput',false);
if matRad_cfg.isMatlab
    nameList = cellfun(@(mc) mc.PropertyList(strcmp({mc.PropertyList.Name}, 'name')).DefaultValue,availableDoseEngines,'UniformOutput',false);
else
    nameList = cellfun(@(mc) mc.PropertyList{strcmp({mc.PropertyList.Name}, 'name')}.DefaultValue,availableDoseEngines,'UniformOutput',false);
end
handleList = cellfun(@(namestr) str2func(namestr),classList,'UniformOutput',false);

end


