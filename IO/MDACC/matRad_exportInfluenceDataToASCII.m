function [ flagSuccess ] = matRad_exportInfluenceDataToASCII(cst,stf,pln,dij )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ExportPath       = '/Users/Hans-PeterWieser/Documents/Heidelberg/MDAnderson2016';
formatSpecIJ     = '%d.%d.%d.%d.%d %.8f\n';  % format in files is:   119.93.207.0.1000 4.064967e-01
formatSpecStruct = '%d,%d,%d\n';
ixExport   = find(~cellfun(@isempty, dij.physicalDose))'; 

rootFolder   = [ExportPath filesep 'matRadExport_' datestr(datetime('now'),'mm-dd-yyyy-HH-MM-SS')];
structFolder = [rootFolder filesep 'structs'];
mkdir(rootFolder);
mkdir(structFolder);

%% export structures
for i = 1:size(cst,1)
      [x,y,z] = ind2sub(dij.dimensions,cst{i,4}{1});
      data = [x,y,z]';
      fileID  = fopen([structFolder filesep lower(cst{i,2}) '.txt'],'w');
      fprintf(fileID,formatSpecStruct,data);
      fclose(fileID);
end


%% export influence data
CntBeamlet = 1;
ixLow      = 1;
for i = 1:numel(ixExport)
   
   
   subFolder = ['scenario_' num2str(i,'%d')];
   mkdir([rootFolder filesep subFolder])
   fullPath = [rootFolder filesep subFolder];
   
   [linIx,beamletIx,dose]   = find(dij.physicalDose{ixExport(i)});
   [x,y,z]                  = ind2sub(dij.dimensions,linIx);
   if isfield(dij,'mLETDose')
      [~,~,LETdose]         = find(dij.mLETDose{ixExport(i)});
   end
   
   for beamIx = 1:pln.numOfBeams
      
      UpperLimit = CntBeamlet + stf(beamIx).totalNumOfBixels - 1;
      ixUpp      = find(beamletIx<=UpperLimit,1,'last');
      vBeamIx    = ones(ixUpp-ixLow+1,1) * beamIx-1;
     
      
      data = [x(ixLow:ixUpp) y(ixLow:ixUpp) z(ixLow:ixUpp) vBeamIx [beamletIx(ixLow:ixUpp)-1] dose(ixLow:ixUpp)]';
      
      fileID     = fopen([fullPath filesep 'Dij_' num2str(beamIx,'%d') '.ASCII'],'w');
      fprintf(fileID,formatSpecIJ,data);
      fclose(fileID);
      
      
      if isfield(dij,'mLETDose')
         data = [x(ixLow:ixUpp) y(ixLow:ixUpp) z(ixLow:ixUpp) vBeamIx beamletIx(ixLow:ixUpp)-1 LETdose(ixLow:ixUpp)./dose(ixLow:ixUpp)]';
         
         fileID = fopen([fullPath filesep 'LETij_' num2str(beamIx,'%d') '.ASCII'],'w');
         fprintf(fileID,formatSpecIJ,data);
         fclose(fileID);
      end
      
      CntBeamlet = CntBeamlet + stf(beamIx).totalNumOfBixels;
      ixLow = ixLow + (ixUpp-ixLow+1);
      
   end
end


flagSuccess = true;
end


