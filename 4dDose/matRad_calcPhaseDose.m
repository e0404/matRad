function [resultGUI, delivery] = matRad_calcPhaseDose();
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
% call
%   
%
% input
%       
%  
% output
%
% comment:
% 
% References
%
% Silke Ulrich Aug 2016
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get data from workspace
multScen       = evalin('base','multScen');
dij            = evalin('base','dij');
stf       = evalin('base','stf');

NumberOfBeams = length(stf);
for i=1:NumberOfBeams
    allocate beam spots to time points
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
        while(v<delivery(i).IES(IES).voxel && j < length(delivery(i).es))  
            while(j < length(delivery(i).es) && delivery(i).np(j)+1 == delivery(i).np(j+1))
                delivery(i).energy(j) = str2num(xPln(i).PTTxPlanMd5.PTTxPlan.Beam.IES{1,IES}.Attributes.energy);
                delivery(i).xpos(j) = str2num(xPln(i).PTTxPlanMd5.PTTxPlan.Beam.IES{1,IES}.Voxel{1,v}.Attributes.x);
                delivery(i).ypos(j) = str2num(xPln(i).PTTxPlanMd5.PTTxPlan.Beam.IES{1,IES}.Voxel{1,v}.Attributes.y);
                delivery(i).particles(j) = str2num(xPln(i).PTTxPlanMd5.PTTxPlan.Beam.IES{1,IES}.Voxel{1,v}.Attributes.particles);
                j = j+1;
                v=v+1;
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

    % find correct j-value for dose calculation with pre-calculated dij
    % Matrix
    % von Lucas export HITXMLPlan : HITXML, x-y-z(beam)
    %                               matRad, x-y(beam)-z -> z-x-y(beam)
    %      Koordinaten              voxel_x = rayPos_bev(3);
    %                               voxel_y = rayPos_bev(1);
         
    delivery(i).stf_rayNum = zeros(length(delivery(i).es), 1); 
    delivery(i).stf_energyIx = zeros(length(delivery(i).es), 1); 
    delivery(i).j = zeros(length(delivery(i).es), 1); 
    for c=1:length(delivery(i).time)-1
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

NumOfPhases = multScen.numOfCtScen;


