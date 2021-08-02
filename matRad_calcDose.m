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

if isfield(pln,'propDoseCalc') 
    
    %get all available engines for given pln struct
    [nameList, ~, handleList] = DoseEngines.matRad_DoseEngine.getAvailableEngines(pln);
    %check if there is an engine field in the pln struct
    if (isfield(pln.propDoseCalc,'engine'))
        
        %check if there is already an object or only a struct 
        if~isstruct(pln.propDoseCalc.engine)        
                     
            %check if it is a dose calc engine and not some other oject
            if isa(pln.propDoseCalc.engine, 'DoseEngines.matRad_DoseEngine')
                %call the calcDose funktion 
                dij = pln.propDoseCalc.engine.calcDose(ct,stf,pln,cst);
            else
                matRad_cfg.dispError('Pln struct contains non valid dose engine: %s', pln.propDoseCalc.engine.name);
            end
                   
        else
            %engine field is a struct
            if (isfield(pln.propDoseCalc.engine ,'name') && any(strcmpi(nameList,pln.propDoseCalc.engine.name)))
                
                engineHandle = handleList{strcmpi(nameList,pln.propDoseCalc.engine.name)};
                engine = engineHandle(ct,stf,pln,cst);
                fields = fieldnames(pln.propDoseCalc.engine);
                for i = 1:length(fields)
                    
                    try
                        engine.(fields{i}) = pln.propDoseCalc.engine.(fields{i});
                        
                    catch ME
                        switch Me.identifier
                            case 'MATALB:UndefinedProperty' 
                                matRad_cfg.dispWarning('ME field:%s %s',fields{i},ME.identifier);
                            otherwise
                                rethrow(ME);
                        end
                    end         
                    
                end
                
                dij = engine.calcDose(ct,stf,pln,cst);
                
                                
            else
                matRad_cfg.dispError('No valid engine name given in pln.propDoseCalc struct');
            end
        
        end
               
        
    else
        %if no engine informations are given use the std engine from matRad_Config 
        
    end

else
    
    matRad_rc.dispError('No propDoseCalc field inside given pln struc.');

end

