function matRad_unitTestTextManipulation()

fid=fopen('matRad_calcPhotonDose.m');
fo=fopen('tempFile1.m','w');
tline = fgetl(fid);

while ischar(tline)
    
    if strfind(tline, 'lateralCutoff = 50')
        fprintf(fo, '%s\n', 'lateralCutoff = 20;');
    else
        fprintf(fo, '%s\n', tline);
    end
    tline = fgetl(fid);
end

fclose(fid);
fclose(fo);

movefile('tempFile1.m', '../matRad_calcPhotonDose.m', 'f');


fid=fopen('matRad_calcParticleDose.m');
fo=fopen('tempFile2.m','w');
tline = fgetl(fid);
while ischar(tline)
    
    if strfind(tline, 'cutOffLevel          = 0.99')
        fprintf(fo, '%s\n', '       cutOffLevel          = 0.8;');
    else
        fprintf(fo, '%s\n', tline);
    end
    tline = fgetl(fid);
end

fclose(fid);
fclose(fo);

movefile('tempFile2.m', '../matRad_calcParticleDose.m', 'f');

fid=fopen('matRad_ipoptOptions.m');
fo=fopen('tempFile3.m','w');
tline = fgetl(fid);

while ischar(tline)
    
    if strfind(tline, 'options.ipopt.max_iter')
        fprintf(fo, '%s\n', 'options.ipopt.max_iter = 10;');
    else
        fprintf(fo, '%s\n', tline);
    end
    tline = fgetl(fid);
end

fclose(fid);
fclose(fo);

movefile('tempFile3.m', '../optimization/matRad_ipoptOptions.m', 'f');
end

