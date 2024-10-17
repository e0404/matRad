function dataOut = matRad_readCsvData(csvFile,cubeDim)
% matRad read TOPAS csv data
%
% call
%   dataOut = matRad_readCsvData(csvFile,cubeDim)
%
% input
%   csvFile: TOPAS csv scoring file
%   cubeDim: size of cube
%
% output
%   dataOut: cube of size cubeDim containing scored values
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2024 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read in csv file as table
dataTable = readtable(csvFile,'ReadVariableNames',false);

% this is the number of ReportQuantities contained in that file
numOfReportQuantities = size(dataTable,2)-3;

% Save all scored quantities as cell array and reshape to cubeDim
dataOut = cell(1,numOfReportQuantities);

%Get the indices and do i,j swap
ixi = dataTable.Var2+1;
ixj = dataTable.Var1+1;
ixk = dataTable.Var3+1;
ix = sub2ind(cubeDim,ixi,ixj,ixk);

for i = 1:numOfReportQuantities
    dataOut{i} = zeros(cubeDim);
    dataOut{i}(ix) = dataTable.(['Var' num2str(i+3)]);
end

end