function [resultGUI, delivery] = matRad_calcPhaseDose(resultGUI, dij, delivery)
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

MOTION = 'linear'; %'CALYPSO' 'ANZAI' 

add = 0;
for i=1:length(delivery)
    delivery(i).phase = zeros(length(delivery(i).es), 1); 
    delivery(i).w = zeros(length(delivery(i).es), 1); 
   %delivery(i).index_j = zeros(length(delivery(i).es), 1); 
    
    % ordne jeden spot sein weight zu
    % [j, index_j] = sort(delivery(i).j);
    % delivery(i).index_j = index_j;
    for l=1:length(delivery(i).j)
        if(isnan(delivery(i).j(l)))
        delivery(i).w(l) = NaN;
    else
        delivery(i).w(l) = resultGUI.w(delivery(i).j(l));  
    end
    end
    %add = max(j);


%%%%%%%%%Breathing motion   
    %Motion assumption1: linear
    if(strcmp(MOTION,'linear'))
        NumOfPhases = size(dij.physicalDose);
        NumOfPhases =  NumOfPhases(1);  %immer richtig?
    offset = delivery(1).offset; %0;
    phaseTime = delivery(1).motionperiod/NumOfPhases; %5/NumOfPhases;
    t=1;
    p = ceil((offset+0.001)/phaseTime);  
    while p > NumOfPhases
        p = p-NumOfPhases;
    end
    
    Time = mod(offset, phaseTime); %offset; 
    if(Time == 0)
        Time = Time + phaseTime;
    end
    while(t < length(delivery(i).time))       
        while(t < length(delivery(i).time) && delivery(i).time(t)/1000000 < Time)  
            delivery(i).phase(t) = p;
            t= t+1;
        end
        p = p+1;
        if(p > NumOfPhases)
            p=1;
        end  
         Time = Time + phaseTime;
    end
      
  
    %Motion assumption3: from Calypso data
    %amplitude based
    % evtl. noch phase base, relative phased based according to Richter?
    elseif(strcmp(MOTION,'CALYPSO'))
    fid = fopen('D:\Matrad\data\4DCT\TKUH005_linux\Motion_005\RawDataTxt\Patient_005\CalypsoLung_005_F01_TrackTarget.txt');
    dataArray = textscan(fid, '%f %f %f %f %d', 'headerlines',1);
    
    motion.time = dataArray{1};
    motion.x = dataArray{2};
    motion.y = dataArray{3};
    motion.z = dataArray{4};
    fclose(fid);   

      
  %Motion assumption3: from Anzai data
    elseif(strcmp(MOTION, 'ANZAI'))
    fid = fopen('D:\Silke\data\CalypsoLunge\Anzai-Daten TKUH006-bisTKUH008\Anzai-Daten TKUH006-bisTKUH008\Anzai_006\Anzai_006_11.txt');
    dataArray = textscan(fid, '%f %f %d', 'headerlines',1); 
    
    motion.time = dataArray{1};
    motion.amplitude = dataArray{2};
   
    figure
    plot (motion.time, motion.amplitude)
    title('Anzai motion')
    xlabel('time [s]')
    ylabel('amplitude')
    
    %todo: Zuordnung welche Phase zu welcher Amplitude
    % amplitude or phase based???
 
end
    %%%%%%%%%%%%%%%END MOTION
    
end %beams



%berechne Dosis für jede Phase  
disp('calc dose in each CT phase')
for p=1:NumOfPhases
    w=zeros(dij.totalNumOfBixels, 1); 
    for i=1:length(delivery)
        t =  find(delivery(i).phase == p  & ~isnan(delivery(i).w));
        w(delivery(i).j(t)) = delivery(i).w(t);       
    end
 
    resultGUI.phaseDose{p} = reshape(dij.physicalDose{p} * w, dij.dimensions);
    
    if isequal(resultGUI.bioParam.type,'MCN_RBExD')
    
        resultGUI.phaseAlphaDose{p} = reshape(dij.mAlphaDose{p} * w, dij.dimensions);
        resultGUI.phaseSqrtBetaDose{p} = reshape(dij.mSqrtBetaDose{p} * w, dij.dimensions);
    
    end
end
  


