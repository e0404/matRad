function dij = matRad_calcDose(ct,stf,pln,cst)
    % matRad dose calculation automaticly creating the appropriate dose engine
    % for the given pln struct and called the associated dose calculation funtion 
    % 
    % call
    %   dij =  matRad_calcDose(ct,stf,pln,cst)
    %   dij =  matRad_calcDose(ct,stf,pln,cst,calcDoseDirect)
    %
    % input
    %   ct:         ct cube
    %   stf:        matRad steering information struct
    %   pln:        matRad plan meta information struct
    %   cst:        matRad cst struct
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


    %get all available engines for given pln struct, could be done conditional
    [nameList, ~, handleList] = DoseEngines.matRad_DoseEngine.getAvailableEngines(pln);

    if(isfield(pln,'propDoseCalc') && isfield(pln.propDoseCalc,'engine')) 


        % check if there is an engine field in the pln struct        
        % check if there is already an object or only a struct 
        if~isstruct(pln.propDoseCalc.engine)        

            % check if it is a dose calc engine and not some other oject
            if isa(pln.propDoseCalc.engine, 'DoseEngines.matRad_DoseEngine')
                % set engine
                engine = pln.propDoseCalc.engine;
            else
                % cancel if the given engine is not a valid dose calc engine
                matRad_cfg.dispError('pln struct contains non valid dose engine: %s! ', pln.propDoseCalc.engine.name);
            end

        else
           % engine field is a struct
            if (isfield(pln.propDoseCalc.engine ,'name') && any(strcmpi(nameList,pln.propDoseCalc.engine.name)))

                engineHandle = handleList{strcmpi(nameList,pln.propDoseCalc.engine.name)};
                engine = engineHandle(ct,stf,pln,cst);
                fields = fieldnames(pln.propDoseCalc.engine);
                
                % name field is no longer needed and would throw an exception
                fields(strcmp(fields, 'name')) = [];
                
                % itterate over all fieldnames and try to set the
                % corresponding properties inside the engine
                for i = 1:length(fields)                   
                    try
                        engine.(fields{i}) = pln.propDoseCalc.engine.(fields{i});
                    
                    % catch exceptions when the engine has no properties,
                    % which are defined in the struct. 
                    % When defining an engine with custom setter and getter
                    % methods, custom exceptions can be caught here. Be
                    % careful with Octave exceptions!
                    catch ME
                        switch ME.identifier
                            case 'MATLAB:noPublicFieldForClass' 
                                matRad_cfg.dispWarning('Problem with given engine struct: %s',ME.message);
                            otherwise
                                matRad_cfg.dispWarning('Problem while setting up engine from struct:%s %s',fields{i},ME.message);
                        end
                    end         

                end


            else
                % if the given engine isn't valid set boolean to initiliaze
                % it at the end of this function
                initDefaultEngine = true;
                matRad_cfg.dispWarning('No valid engine name for given radiation mode given in pln.propDoseCalc struct. Trying to use default dose calc engine from MatRad_Config.');
                
            end
        end

    else
        % if the given engine isn't valid set boolean to initiliaze
        % it at the end of this function
        initDefaultEngine = true;
        matRad_cfg.dispWarning('No engine or propDoseCalc field inside given pln struct. Trying to use default dose calc engine from MatRad_Config.');

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
            engine = engineHandle(ct,stf,pln,cst);         
            matRad_cfg.dispWarning('Using %s Engine as default engine!', engine.name);

        else
            
            matRad_cfg.dispError('No valid default engine found for radiation mode.');   
        
        end
   
        
    end
    %call the calcDose funktion 
    dij = engine.calcDose(ct,stf,pln,cst);

end
