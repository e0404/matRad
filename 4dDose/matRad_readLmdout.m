function delivery = matRad_readLmdout(dij, stf, FileName, FileName_lmdout)
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
%   delivery 

% comment:
% dij and stf are needed in workspace
% At the moment a lot of variables are written in delivery struct, some of
% them might be deleted 
%   necessary: time ,xpos, ypos, energy, j


% References
%
%
% Silke Ulrich Aug 2016
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%disp('read in delivery information')

if(nargin <4)
    FileName_lmdout = FileName;
end

delivery = struct;

NumberOfBeams = length(stf);

for i=1:NumberOfBeams
    n = num2str(i-1);
    PlnFile = ['PBP_0' n '_' FileName '.xml'];
    LmdoutFile = ['D_0' n '_' FileName_lmdout '.lmdout'];
    
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
    formatSpec = '%f%d%d%d%d%d%d%d%d%*s%*s%[^\n\r]';

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
     
    %allocate beam spots to time points
    delivery(i).energy = zeros(length(delivery(i).es), 1); 
    delivery(i).xpos = zeros(length(delivery(i).es), 1);
    delivery(i).ypos = zeros(length(delivery(i).es), 1);
    delivery(i).particles = zeros(length(delivery(i).es), 1);  % löschen???
    j=1;
    for IES = 1:delivery(i).energies
        delivery(i).IES(IES).energy = str2num(xPln(i).PTTxPlanMd5.PTTxPlan.Beam.IES{1,IES}.Attributes.energy);
         v=1;
         jold = j;
        %while(delivery(i).eop(j) == delivery(i).eop(jold))  % geht nicht da eop auch hochgezählt wird wenn Energie nicht geändert wird ???
        while(v<=delivery(i).IES(IES).voxel && j < length(delivery(i).es))  
            if(delivery(i).IES(IES).voxel == 1) %dann andere Struktur in xPln
                delivery(i).energy(j) = str2num(xPln(i).PTTxPlanMd5.PTTxPlan.Beam.IES{1,IES}.Attributes.energy);
                delivery(i).xpos(j) = str2num(xPln(i).PTTxPlanMd5.PTTxPlan.Beam.IES{1,IES}.Voxel.Attributes.x);
                delivery(i).ypos(j) = str2num(xPln(i).PTTxPlanMd5.PTTxPlan.Beam.IES{1,IES}.Voxel.Attributes.y);
                delivery(i).particles(j) = str2num(xPln(i).PTTxPlanMd5.PTTxPlan.Beam.IES{1,IES}.Voxel.Attributes.particles);
                j = j+1;
                v=v+1;
            else
                while(j < length(delivery(i).es) && delivery(i).np(j)+1 == delivery(i).np(j+1))
                    delivery(i).energy(j) = str2num(xPln(i).PTTxPlanMd5.PTTxPlan.Beam.IES{1,IES}.Attributes.energy);
                    delivery(i).xpos(j) = str2num(xPln(i).PTTxPlanMd5.PTTxPlan.Beam.IES{1,IES}.Voxel{1,v}.Attributes.x);
                    delivery(i).ypos(j) = str2num(xPln(i).PTTxPlanMd5.PTTxPlan.Beam.IES{1,IES}.Voxel{1,v}.Attributes.y);
                    delivery(i).particles(j) = str2num(xPln(i).PTTxPlanMd5.PTTxPlan.Beam.IES{1,IES}.Voxel{1,v}.Attributes.particles);
                    j = j+1;
                    v=v+1;
                end
            end
            while(j < length(delivery(i).es) && delivery(i).np(j) == delivery(i).np(j+1))
                delivery(i).energy(j) = NaN;
                delivery(i).xpos(j) =  NaN;
                delivery(i).ypos(j) =  NaN;
                delivery(i).particles(j) =  NaN;
                j = j+1;
            end
        end  
    end %IES
    
    delivery(i).energy(j) = NaN;
    delivery(i).xpos(j) =  NaN;
    delivery(i).ypos(j) =  NaN;
    delivery(i).particles(j) =  NaN;

    % find correct j-value for dose calculation with pre-calculated dij
    % Matrix
    % von Lucas export HITXMLPlan : HITXML, x-y-z(beam)
    %                               matRad, x-y(beam)-z -> z-x-y(beam)
    %      Koordinaten              voxel_x = rayPos_bev(3);
    %                               voxel_y = rayPos_bev(1);
         
    delivery(i).stf_rayNum = zeros(length(delivery(i).es), 1); 
    delivery(i).stf_energyIx = zeros(length(delivery(i).es), 1); 
    delivery(i).j = zeros(length(delivery(i).es), 1); 
    for c=1:length(delivery(i).time)
        if(~isnan(delivery(i).xpos(c)))
            delivery_pos = [delivery(i).xpos(c); delivery(i).ypos(c)]*ones(1,stf(i).numOfRays);
        
            stf_pos =reshape([stf(i).ray.rayPos_bev],3,[]);
            rayNum = find(sum(stf_pos([3 1],:) == delivery_pos) == 2);
        
            energyIx = find(abs(stf(i).ray(rayNum).energy-delivery(i).energy(c))<0.1);

            j = find(dij.rayNum == rayNum & dij.bixelNum == energyIx & dij.beamNum == i);
            
            delivery(i).stf_rayNum(c) = rayNum;
            delivery(i).stf_energyIx(c) = energyIx;
            delivery(i).j(c) = j;
        else
            delivery(i).stf_rayNum(c) = NaN;
            delivery(i).stf_energyIx(c) = NaN;
            delivery(i).j(c) = NaN;
        end
    end    
    
end %beams












