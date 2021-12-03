function data = readCsvData(csvFile,cubeDim)
    data = zeros(cubeDim(2),cubeDim(1),cubeDim(3));
    fID = fopen(csvFile,'r');
    dataCsv = textscan(fID,'%d %d %d %f','Delimiter',',','CommentStyle','#','CollectOutput',true);
    fclose(fID);
    ix = sub2ind([cubeDim(1) cubeDim(2) cubeDim(3)],dataCsv{1}(:,2)+1,dataCsv{1}(:,1)+1,dataCsv{1}(:,3)+1);
    data(ix) = dataCsv{2};
end