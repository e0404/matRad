function data = readBinData(binFile,cubeDim)
    fID = fopen(binFile);
    data = fread(fID,prod(cubeDim),'double');
    fclose(fID);
    data = reshape(data,cubeDim(2),cubeDim(1),cubeDim(3));
    data = permute(data,[2 1 3]);
end