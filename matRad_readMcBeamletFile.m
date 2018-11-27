function matRad_readMcBeamletFile(filename)

fileHandle = fopen(filename,'r');

numOfBeams = fread(fileHandle,1,'int');

for i = 1:numOfBeams
    currSourcePos = fread(fileHandle,3,'double');
    numOfBeamlets = fread(fileHandle,1,'int');
    for j = 1:numOfBeamlets
        currRayCorner1 = fread(fileHandle,3,'double');
        currRayCorner2 = fread(fileHandle,3,'double');
        currRayCorner3 = fread(fileHandle,3,'double');
        currRayCorner4 = fread(fileHandle,3,'double')
    end
end

fclose(fileHandle);