function valueResult = matRad_readDataFromText(fileName, textBeforeData_1,textBeforeData_2, textBehindData, numberOfColumns)

%% Function to read data between two known lines of text
% Input:    fileName
%           textBeforeData_1  -   last line or text fragment before data
%           textBeforeData_2  -   specific line/fragment before data (if 
%                                 not applicable set to char(''))
%           textBehindData    -   first line behind data
%           numberOfColumns   -   optional reordering into given # of 
%                                 colums    
%
% Author: Lucas Sommer, 10.05.2020 (lucas.sommer@tum.de)

%% Open text file
fid_read = fopen(fileName, 'rt');

%% Find last line before data
lineDummy = fgetl(fid_read);
while(~strcmp(lineDummy(1,1:min(length(lineDummy),length(textBeforeData_1))), textBeforeData_1(1,:)))
    clear lineDummy
    lineDummy = fgetl(fid_read);
end

if ~isempty(textBeforeData_2)
    while(~strcmp(lineDummy(1,1:min(length(lineDummy),length(textBeforeData_2))), textBeforeData_2(1,:)))
        clear lineDummy
        lineDummy = fgetl(fid_read);
    end
end

%% Read data
% valueResult = char(''); %fgetl(fid_read);
% dummyTestEndData = fgetl(fid_read); %valueResult;
% while(~strcmp(dummyTestEndData(1,1:min(length(dummyTestEndData), length(textBehindData))), textBehindData(1,:)))
%     if ~feof(fid_read)
%     valueResult = [valueResult ' ' dummyTestEndData];
%     dummyTestEndData = fgetl(fid_read);
%     elseif feof(fid_read)
%         valueResult = [valueResult ' ' dummyTestEndData];
%         break
%     end
% end
% 
% valueResult = sscanf(valueResult, '%f');

valueResult_char = char('');
valueResult_sparse = spalloc(512^3,1,1);
sparceCounter = 1;

dummyTestEndData = fgetl(fid_read);
while(~strcmp(dummyTestEndData(1,1:min(length(dummyTestEndData), length(textBehindData))), textBehindData(1,:)))
    if ~feof(fid_read)
        valueResult_char = dummyTestEndData;
        valueDummy = sscanf(valueResult_char, '%f');
        valueResult_sparse(sparceCounter:sparceCounter+length(valueDummy)-1,1) = valueDummy(:);
        sparceCounter = sparceCounter+length(valueDummy);       
        dummyTestEndData = fgetl(fid_read);
    elseif feof(fid_read)
        valueResult_char = dummyTestEndData;
        valueDummy = sscanf(valueResult_char, '%f');
        valueResult_sparse(sparceCounter:sparceCounter+length(valueDummy)-1,1) = valueDummy(:);
        sparceCounter = sparceCounter+length(valueDummy);       
        break
    end
end

valueResult = valueResult_sparse(1:sparceCounter-1,1);

%% Close text file
fclose(fid_read);

%% Reorder (optional)
if (nargin > 4) && (numberOfColumns <= 3)
    switch numberOfColumns
        case 1
            disp('Data output given in one column.')
        case 2
            if ~mod(length(valueResult_sparse),2)
                valueResult_sparse = reshape(valueResult_sparse, 2, length(valueResult_sparse)/2);
                disp('Data output given in two columns.')
            elseif mod(length(valueResult_sparse),2)
                disp('Number of values does not match wished number of columns, reordering not possible!')
            end
        case 3
            if ~mod(length(valueResult_sparse),3)
                valueResult_sparse = reshape(valueResult_sparse, 3, length(valueResult_sparse)/3);
                disp('Data output given in three column.')
            elseif mod(length(valueResult_sparse),3)
                disp('Number of values does not match wished number of columns, reordering not possible!')
            end
    end
elseif(nargin > 4) && (numberOfColumns > 3)
    disp('Too many columns were selected, reordering not possible!')
end

end