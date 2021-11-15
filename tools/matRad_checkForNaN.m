function outArray = matRad_checkForNaN(structIn)
fn = fieldnames(structIn);
counter = 1;
for k = 1:numel(fn)
    if isnumeric(structIn.(fn{k}))
        outArray{counter,1} = fn{k};
        outArray{counter,2} = any(isnan(structIn.(fn{k})(:)));
        counter = counter + 1;
    end
end

end