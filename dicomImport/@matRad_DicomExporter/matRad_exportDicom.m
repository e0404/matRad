function obj = matRad_exportDicom(obj)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%Name of Patient & Study
%CT Series
obj = matRad_exportDicomCt(obj);

%RTStruct Series

if ~isempty(obj.cst)
   obj = matRad_exportDicomRTStruct(obj);    
end

if ~isempty(obj.pln)
    obj = matRad_exportDicomRTPlan(obj);
end

if ~isempty(obj.resultGUI)
    obj = matRad_exportDicomRTDoses(obj);
end




%Pln Series

%Dose Series
end

