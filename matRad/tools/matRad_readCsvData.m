function dataOut = matRad_readCsvData(csvFile,cubeDim)

% Instance of MatRad_Config class
matRad_cfg = MatRad_Config.instance();

% Read in csv file as table
dataTable = readtable(csvFile,'ReadVariableNames',false);

% Check if consistent with cubeDim
if rem(size(dataTable,1),prod(cubeDim))==0
    % this is the number of ReportQuantities contained in that file
    numOfReportQuantities = size(dataTable,2)-3;

    % Save all scored quantities as cell array and reshape to cubeDim
    dataOut = cell(1,numOfReportQuantities);
    for i = 1:numOfReportQuantities
        dataOut{i} = reshape(dataTable.(['Var' num2str(i+3)]),cubeDim(2),cubeDim(1),cubeDim(3));
        dataOut{i} = permute(dataOut{i},[2 1 3]);
    end
else
    matRad_cfg.dispError('bin data contains an odd number of entries.')
end

end