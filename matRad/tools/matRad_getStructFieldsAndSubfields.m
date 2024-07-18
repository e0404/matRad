function fieldStruct = matRad_getStructFieldsAndSubfields(s)

    fieldStruct = fieldnames(s);
    currLevelFields  = fieldStruct;
    nFields  = numel(fieldStruct);

    for subIdx=1:nFields
        currField = s.(currLevelFields{subIdx});
        if isstruct(currField)
            subFields = matRad_getStructFieldsAndSubfields(currField);
            fieldStruct = [fieldStruct; cellfun(@(sF) [fieldStruct{subIdx}, '.', sF], subFields, 'UniformOutput',false)];
        end
    end
end