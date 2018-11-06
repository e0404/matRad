function matRad_unitTestTextManipulation(filename, string1, string2)

fid=fopen(filename);
fo=fopen('tempFile.m','w');
tline = fgetl(fid);

while ischar(tline)
    
    if strfind(tline, string1)
        fprintf(fo, '%s\n', string2);
    else
        fprintf(fo, '%s\n', tline);
    end
    tline = fgetl(fid);
end

fclose(fid);
fclose(fo);

movefile('tempFile.m', '../matRad_calcPhotonDose.m', 'f');
end

