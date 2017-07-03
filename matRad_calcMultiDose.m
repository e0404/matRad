function [DoseMat] = matRad_calcMultiDose(TracMat,baseData,sig)


TracMatc = TracMat;
cutt = any(TracMat,1);
TracMatc( :, ~cutt) = [];
cutt2 = any(TracMatc,2);
TracMatc(~cutt2,:) = [];

dimTM = size(TracMatc);

[x,~]=meshgrid(1:dimTM(2),1:dimTM(1));
[x1,~]=meshgrid(1:dimTM(2),1:length(baseData.Z));

interpData = griddata(repmat(baseData.depths,[1 dimTM(2)]),x1,...
    repmat(baseData.Z,[1 dimTM(2)]),TracMatc,x);

TracMatf = zeros(size(TracMat));
TracMatf(cutt2,cutt) = interpData;


