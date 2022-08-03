function [nameList, classList, handleList] = getAvailableEngines(pln,optionalPath)
    % Returns a list of names and coresponding handle for available dose calc engines
    %   Returns all dose calc engines in the package when no arg is
    %   given. If no engines are found return gonna be empty.
    %
    % call:
    %   [nameList, handleList] = DoseEngines.matRad_DoseEngine.getAvailableEngines(pln,optional_path)
    %
    % input:
    %   pln: containing proposed dose calc and machine file informations
    %   optionalPath: searches for dose calc engines in given
    %
    % returns:
    %   nameList: cell-array conatining readable names for engines
    %   classList: cell-array conatining full classnamens for available engines
    %   handleList: cell-array containing function-handles to
    %                  available engines constructor (call the included handle by adding Parentheses e.g. handleList{1}())

    matRad_cfg = MatRad_Config.instance();


    nameList = {};
    classList = {};
    handleList = {};

    switch nargin

        case 0
            % meta.package works perfectly in matlab but won't work in Octave
            % because Octave doesn't fully support package folder and class definition,
            % so we use a somewhat unclean work around in Octave which propably won't work withg class folder.
            % So when using engine classes in octave it's better to initionlize them from hand
            if matRad_cfg.isOctave
                [nameList, classList, handleList] = matRad_getAvailableEnginesOctave();
            else
                % get all meta classes in the DoseEngines package
                mp = meta.package.fromName('DoseEngines');
                mc_list = mp.ClassList;

                % itterate through the meta classes in the package
                for i = 1:length(mc_list)

                    if matRad_cfg.isOctave
                        [~,className] = fileparts(files(i).name);
                        % get meta class with package name infront
                        mc = meta.class.fromName(['DoseEngines.' className]);
                    else
                        mc = mc_list(i);
                    end
                    % check for the isCalcEngine property,
                    % could be done cleaner with the superclasses method
                    % which sadly isn't available in octave
                    [~,mc_isEngineIdx] = ismember('isDoseEngine', {mc.PropertyList.Name});
                    if ((~isempty(mc.SuperclassList) && any(strcmp(mc.SuperclassList.Name, 'DoseEngines.matRad_DoseEngine'))) || ...
                            (mc_isEngineIdx && mc.PropertyList(mc_isEngineIdx).DefaultValue))

                        handleList{end+1} = str2func(mc.Name);
                        classList{end+1} = mc.Name;

                        %get readable name from metaclass properties
                        nameProp = mc.PropertyList(strcmp({mc.PropertyList.Name}, 'name'));
                        if (nameProp.Abstract)
                            % ND -> not defined meaning abstract class
                            % without a name
                            nameList{end+1} = 'ND';
                        else
                            nameList{end+1} = nameProp.DefaultValue;
                        end


                    end
                end

            end
        case 1
            if matRad_cfg.isOctave
                [nameList, classList, handleList] = matRad_getAvailableEnginesOctave(pln);
            else
                mp = meta.package.fromName('DoseEngines');
                mc_list = mp.ClassList;
                % itterate through the meta classes in the package
                for i = 1:length(mc_list)

                    mc = mc_list(i);
                    % skip class if abstract
                    if~(mc.Abstract)
                        % check for the isCalcEngine property,
                        % could be done cleaner with the superclasses method
                        % which sadly isn't available in octave
                        [~,mc_isEngineIdx] = ismember('isDoseEngine', {mc.PropertyList.Name});
                        if ((~isempty(mc.SuperclassList) && any(strcmp(mc.SuperclassList.Name, 'DoseEngines.matRad_DoseEngine'))) || ...
                                (mc_isEngineIdx && mc.PropertyList(mc_isEngineIdx).DefaultValue))

                            % get radiation mode from meta class property
                            [~, loc] = ismember('possibleRadiationModes', {mc.PropertyList.Name});
                            propValue = mc.PropertyList(loc).DefaultValue;

                            if(any(strcmp(propValue, pln.radiationMode)))
                                % get radiation mode from the in pln proposed basedata machine file
                                machineMode = DoseEngines.matRad_DoseEngine.loadMachine(pln).meta.radiationMode;

                                % add current class to return lists if the
                                % radiation mode is compatible
                                if(any(strcmp(propValue, machineMode)))
                                    handleList{end+1} = str2func(mc.Name);
                                    classList{end+1} = mc.Name;

                                    %get readable name from metaclass properties
                                    nameProp = mc.PropertyList(strcmp({mc.PropertyList.Name}, 'name'));
                                    if (nameProp.Abstract)
                                        % ND -> not defined meaning abstract class
                                        % without a name
                                        nameList{end+1} = 'ND';
                                    else
                                        nameList{end+1} = nameProp.DefaultValue;
                                    end

                                end

                            end

                        end

                    end
                end

            end

        case 2
            if matRad_cfg.isOctave
                [nameList, classList, handleList] = matRad_getAvailableEnginesOctave(pln,optionalPath);
            else
                % check if path is valid and add it to the current
                % matlab path
                if(isfolder(optionalPath))
                    addpath(optionalPath);
                end

                %get all MATLAB relevant files
                pathContent = what(optionalPath);
                %concatenate all .m files and class folder
                files = vertcat(pathContent.m, pathContent.classes);
                for i = 1:length(files)
                    [~,className] = fileparts(files{i});
                    mc = meta.class.fromName(className);
                    if (~isempty(mc) && ~(mc.Abstract))

                        % get the index of the is engine property
                        [~,mc_isEngineIdx] = ismember('isDoseEngine', {mc.PropertyList.Name});

                        if ((~isempty(mc.SuperclassList) && any(strcmp(mc.SuperclassList.Name, 'DoseEngines.matRad_DoseEngine'))) || ...
                                (mc_isEngineIdx && mc.PropertyList(mc_isEngineIdx).DefaultValue))

                            % get radiation mode from meta class property
                            [~, loc] = ismember('possibleRadiationModes', {mc.PropertyList.Name});
                            propValue = mc.PropertyList(loc).DefaultValue;

                            if(any(strcmp(propValue, pln.radiationMode)))
                                % get radiation mode from the in pln proposed basedata machine file
                                machineMode = DoseEngines.matRad_DoseEngine.loadMachine(pln).meta.radiationMode;

                                % add current class to return lists if the
                                % radiation mode is compatible
                                if(any(strcmp(propValue, machineMode)))
                                    handleList{end+1} = str2func(mc.Name);
                                    classList{end+1} = mc.Name;

                                    %get readable name from metaclass properties
                                    nameProp = mc.PropertyList(strcmp({mc.PropertyList.Name}, 'name'));
                                    if (nameProp.Abstract)
                                        % ND -> not defined meaning abstract class
                                        % without a name
                                        nameList{end+1} = 'ND';
                                    else
                                        nameList{end+1} = nameProp.DefaultValue;
                                    end
                                end

                            end

                        end
                    end
                end
            end
    end
end


