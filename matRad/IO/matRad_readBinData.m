function dataOut = matRad_readBinData(binFile,cubeDim)

% Instance of MatRad_Config class
matRad_cfg = MatRad_Config.instance();

% Read in bin file
fID = fopen(binFile);
data = fread(fID,inf,'double');
fclose(fID);

% Check if consistent with cubeDim
if rem(numel(data),prod(cubeDim))==0
    % this is the number of ReportQuantities contained in that file
    numOfReportQuantities = numel(data)/prod(cubeDim);

    % Save all scored quantities as cell array and reshape to cubeDim
    dataOut = cell(1,numOfReportQuantities);
    for i = 1:numOfReportQuantities
        dataOut{i} = reshape(data(i:numOfReportQuantities:end),cubeDim(2),cubeDim(1),cubeDim(3));
        dataOut{i} = permute(dataOut{i},[2 1 3]);
    end
else
    matRad_cfg.dispError('bin data contains an odd number of entries.')
end

end
