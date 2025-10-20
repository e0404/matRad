function tallyData = matRad_readDataFromText_TMESHvBioOpti(fileName, tallyIdentifier, numberOfColumns)

if strcmpi(tallyIdentifier, 'TMESH3')
    content = fileread(fileName) ;
    
    % Find data for TMESH typ 3 tally
    dataBegin_tally   = strfind( content, 'tally    3' );
    dataBegin = strfind( content, 'vals' )+8;
    dataEnd = strfind(content, 'tfc' )-1;
    
    if ~isempty(dataEnd)        
        if ~(dataBegin(end-1)<dataBegin_tally)||~(dataBegin(end)>dataBegin_tally)
            error('TMESH data location does not match expected location at end of tally file.')
        end
    elseif isempty(dataEnd)
        disp('No MCNP tally for RBE calculation.')
    end
    
    tallyData = sscanf(content(dataBegin(end):end), '%f');
    
    %% Reorder (optional)
    if exist('numberOfColumns')
        switch numberOfColumns
            case 1
                disp('Data output given in one column.')
            case 2
                if ~mod(length(tallyData),2)
                    tallyData = reshape(tallyData, 2, length(tallyData)/2);
                    disp('Data output given in two columns.')
                elseif mod(length(tallyData),2)
                    disp('Number of values does not match wished number of columns, reordering not possible!')
                end
            case 3
                if ~mod(length(tallyData),3)
                    tallyData = reshape(tallyData, 3, length(tallyData)/3);
                    disp('Data output given in three column.')
                elseif mod(length(tallyData),3)
                    disp('Number of values does not match wished number of columns, reordering not possible!')
                end
        end
    end
    
    
elseif strcmpi(tallyIdentifier, 'F6heating4RBEcalc')
    content = fileread(fileName) ;
    
    dataBegin = strfind( content, 'vals' )+8;
    dataBegin = dataBegin(1:end-1); % Attention: last entry in dataBegin belongs to TMESH tally for total dose tallying
    dataEnd = strfind(content, 'tfc' )-1;
    
%     f6TallyList = [1016, 1026, 1036, 1046, 1056, ...
%         2016, 2026 2036, 2046, 2056, ...
%         3016, 3026 3036, 3046, 3056, ...
%         4016, 4026 4036, 4046, 4056, ...
%         5016, 5026 5036, 5046, 5056, ...
%         6016, 6026 6036, 6046, 6056, ...
%         7016, 7026 7036, 7046, 7056, ...
%         8016 8036, 8056, ...
%         9016, 1116];

% List w/o anoxic tallies
    f6TallyList = [1016, 1026, 1036, ...
        2016, 2026 2036, ...
        3016, 3026 3036, ...
        4016, 4026 4036, ...
        5016, 5026 5036, ...
        6016, 6026 6036, ...
        8016 8036];

    for counter=1:size(dataEnd,2)
        dataBegin_tally(counter)   = strfind( content, ['tally ', int2str(f6TallyList(counter))] );
    end
    
    
    [~, indexTallyData] = sort(dataBegin_tally);
    
    for counter =1:size(dataEnd,2)
        tallyData.(['tally', int2str(f6TallyList(indexTallyData((counter))))]) = sscanf(content(dataBegin(counter):dataEnd(counter)), '%f');
        
        %% Reorder (optional)
        if exist('numberOfColumns')
            switch numberOfColumns
                case 1
                    disp('Data output given in one column.')
                case 2
                    if ~mod(length(tallyData.(['tally', int2str(f6TallyList(indexTallyData((counter))))])),2)
                        tallyData.(['tally', int2str(f6TallyList(indexTallyData((counter))))]) = reshape(tallyData.(['tally', int2str(f6TallyList(indexTallyData((counter))))]), 2, length(tallyData.(['tally', int2str(f6TallyList(indexTallyData((counter))))]))/2);
                        disp('Data output given in two columns.')
                    elseif mod(length(tallyData.(['tally', int2str(f6TallyList(indexTallyData((counter))))])),2)
                        disp('Number of values does not match wished number of columns, reordering not possible!')
                    end
                case 3
                    if ~mod(length(tallyData.(['tally', int2str(f6TallyList(indexTallyData((counter))))])),3)
                        tallyData.(['tally', int2str(f6TallyList(indexTallyData((counter))))]) = reshape(tallyData.(['tally', int2str(f6TallyList(indexTallyData((counter))))]), 3, length(tallyData.(['tally', int2str(f6TallyList(indexTallyData((counter))))]))/3);
                        disp('Data output given in three column.')
                    elseif mod(length(tallyData.(['tally', int2str(f6TallyList(indexTallyData((counter))))])),3)
                        disp('Number of values does not match wished number of columns, reordering not possible!')
                    end
            end
        end
        
    end
    
    
end

end