function [data] = loadSource(file)

%% Source
global mm cm m mm2 cm2 m2 mm3 cm3 m3 deg uGy mGy cGy Gy seconds minutes h d U Ci mCi Bq



data = load(file);
uNames = fieldnames(data.Units);

for un =1:length(uNames)
    unit   = 1;
    if ~iscell(data.Units.(uNames{un})) && ~isstruct(data.Units.(uNames{un}))
        indStr = strfind(data.Units.(uNames{un}),'^-1');
        
        if indStr
            actUnit = data.Units.(uNames{un})(1:indStr-1);
        else
            actUnit = data.Units.(uNames{un});
        end
        actUnit = sscanf(actUnit,'%s');
        if exist(actUnit)
            if indStr
                unit = unit/eval(actUnit);
            else
                unit = unit*eval(actUnit);
            end
        else
            error('Unable to find the unit')
        end
        data.(uNames{un}) = data.(uNames{un})*unit;
    elseif iscell(data.Units.(uNames{un}))
        for i=1:length(data.Units.(uNames{un}))
            indStr = findstr(data.Units.(uNames{un}){i},'^-1');
            if indStr
                actUnit = data.Units.(uNames{un}){i}(1:indStr-1);
            else
                actUnit = data.Units.(uNames{un}){i};
            end
            actUnit = sscanf(actUnit,'%s');
            if exist(actUnit)
                if indStr
                    unit = unit/eval(actUnit);
                else
                    unit = unit*eval(actUnit);
                end
            else
                error('Unable to find the unit')
            end
        end
        data.(uNames{un}) = data.(uNames{un})*unit;
    elseif isstruct(data.Units.(uNames{un}))
        sNames = fieldnames(data.Units.(uNames{un}));
           
        for sn =1:length(sNames)
             unit   = 1;
            if ~iscell(data.Units.(uNames{un}).(sNames{sn}))
                indStr = strfind(data.Units.(uNames{un}).(sNames{sn}),'^-1');
                
                if indStr
                    actUnit = data.Units.(uNames{un}).(sNames{sn})(1:indStr-1);
                else
                    actUnit = data.Units.(uNames{un}).(sNames{sn});
                end
                actUnit = sscanf(actUnit,'%s');
                if exist(actUnit)
                    if indStr
                        unit = unit/eval(actUnit);
                    else
                        unit = unit*eval(actUnit);
                    end
                else
                    error('Unable to find the unit')
                end
                data.(uNames{un}).(sNames{sn}) = data.(uNames{un}).(sNames{sn})*unit;
            else
                for i=1:length(data.Units.(uNames{un}).(sNames{sn}))
                    indStr = findstr(data.Units.(uNames{un}).(sNames{sn}){i},'^-1');
                    if indStr
                        actUnit = data.Units.(uNames{un}).(sNames{sn}){i}(1:indStr-1);
                    else
                        actUnit = data.Units.(uNames{un}).(sNames{sn}){i};
                    end
                    actUnit = sscanf(actUnit,'%s');
                    if exist(actUnit)
                        if indStr
                            unit = unit/eval(actUnit);
                        else
                            unit = unit*eval(actUnit);
                        end
                    else
                        error('Unable to find the unit')
                    end
                end
                data.(uNames{un}).(sNames{sn}) = data.(uNames{un}).(sNames{sn})*unit;
            end
        end
        
    end
    
end

