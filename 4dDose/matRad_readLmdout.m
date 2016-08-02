function delivery = matRad_readLmdout(FileName,NumberOfBeams)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reads in xml Plan file and lmdout file
% one file per beam
%   XML file:  PBP_BeamNumber_FileName.xml
%   LMDout:    D_BeamNumber_FileName.lmdout

% call
%   
%
% input
%   FileName:      
%   NumberOfBeams:      
%  
% output

% comment:
% At the moment a lot of variables are written in delivery struct, some of
% them might be deleted 
%   necessary: time ,xpos, ypos, energy, j

% References
%
%
% Silke Ulrich Aug 2016
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delivery = struct;

for i=1:NumberOfBeams
    n = num2str(i-1);
    PlnFile = ['PBP_0' n '_' FileName '.xml'];
    LmdoutFile = ['D_0' n '_' FileName '.lmdout'];
    
    %Read in XML Plan file  (energy, spot position)
    if exist(PlnFile, 'file') == 2
        xPln(i) = xml2struct(PlnFile);

        delivery(i).gantryAngle = xPln(i).PTTxPlanMd5.PTTxPlan.Beam.Gantry.Attributes.angle;
        delivery(i).energies = length(xPln(i).PTTxPlanMd5.PTTxPlan.Beam.IES);
        voxel = 0;
        for IES = 1:delivery(i).energies
            voxel = voxel + length(xPln(i).PTTxPlanMd5.PTTxPlan.Beam.IES{1,IES}.Voxel);
            delivery(i).IES(IES).voxel = length(xPln(i).PTTxPlanMd5.PTTxPlan.Beam.IES{1,IES}.Voxel);
        end
        delivery(i).voxel = voxel;  %beam spots
    
    else
        warningMessage = sprintf('Warning: Plnfile does not exist');
        uiwait(msgbox(warningMessage));
    end
    
    % Read in dose delivery simulation file
    if exist(LmdoutFile, 'file') ==2
        fid = fopen(LmdoutFile,'r');
    
    delimiter = {' '};
    startRow = 11;
    endRow = inf;
    formatSpec = '%d%d%d%d%d%d%d%d%d%*s%*s%[^\n\r]';

    dataArray = textscan(fid, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
    
    delivery(i).time = dataArray{1};
    delivery(i).np = dataArray{3};
    delivery(i).bs = dataArray{4};
    delivery(i).es = dataArray{5};
    delivery(i).eop = dataArray{6};

    fclose(fid);
    else
        warningMessage = sprintf('Warning: Lmdoutfile does not exist');
        uiwait(msgbox(warningMessage));
    end
      
end %beams












