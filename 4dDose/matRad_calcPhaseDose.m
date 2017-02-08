function [resultGUI, delivery] = matRad_calcPhaseDose(delivery)
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

MOTION = 'linear'; %'CALYPSO'

% get data from workspace
multScen       = evalin('base','multScen');
resultGUI      = evalin('base','resultGUI');
dij            = evalin('base','dij');

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
        delivery(i).w(l) = resultGUI.finalw(delivery(i).j(l));  %habe in delivery schon particles -> sollte gleich sein!!!
    end
    end
    %add = max(j);


%%%%%%%%%Breathing motion   
    %Motion assumption1: linear
    if(MOTION == 'linear')
    NumOfPhases = multScen.numOfCtScen;
    %Annahme 5s und linear 
    offset = delivery(1).offset; %0;
    phaseTime = delivery(1).motionperiod/NumOfPhases; %5/NumOfPhases;
    j=1;
    p = 1;
    Time = offset;
    while(j < length(delivery(i).time))
        Time = Time + phaseTime;
        while(j < length(delivery(i).time) && delivery(i).time(j)/1000000 < Time)  
            delivery(i).phase(j) = p;
            j= j+1;
        end
        p = p+1;
        if(p > NumOfPhases)
            p=1;
        end
    end
    %Motion assumption2: periodic
    
    
    %Motion assumption3: from Calypso data
    %amplitude based
    % evtl. noch phase base, relative phased based according to Richter?
    elseif(MOTION == 'CALYPSO')
    fid = fopen('D:\Matrad\data\4DCT\TKUH005_linux\Motion_005\RawDataTxt\Patient_005\CalypsoLung_005_F01_TrackTarget.txt');
    dataArray = textscan(fid, '%f %f %f %f %d', 'headerlines',1);
    
    motion.time = dataArray{1};
    motion.x = dataArray{2};
    motion.y = dataArray{3};
    motion.z = dataArray{4};
    fclose(fid);
    
    %nehme nicht gesamt Aufnahme Calypso - evtl. anpassen für anderes
    %motion file
    motion.time = motion.time(500:4000);
    motion.x = motion.x(500:4000);
    motion.y = motion.y(500:4000);
    motion.z = motion.z(500:4000);
    
    %nehme motion in y -Richtung  - 3D Vektor???
    motion.r = motion.y;
    
    % zur Zeit erste CT Phase entspricht end-exhale Positon, also max
    % motion (ode rmin?????????)
    motion_min = min(motion.r);
    motion_max = max(motion.r);    
    motion_absrange = abs(motion_max) + abs(motion_min);
    motion_phaserange = 2*motion_absrange/NumOfPhases; 
    
%     motion.min(1:size(motion.r)) = abs(motion_min);
%     motion.min = motion.min';    
%     motion.r = motion.r + motion.min;

    imin = motion_max - motion_phaserange/2;     imax = imin;
    motion.phase{1,1}=[imin, imax];
    motion.phase{2,1}=0;
    start = true;
    for p=2:NumOfPhases
        if(~start)
           imin = imax;
           imax = imax + motion_phaserange;
           if(imax > motion_max)
            imax = motion_max -imax;
            start = true;
            end                 
        else
        imin = imax;
        imax = imax - motion_phaserange;
        if(imax < motion_min)
            imax = motion_min + (motion_min - imax); %?
            start = false;
        end
        end
        
        motion.phase{1,p}=[imin, imax];
        if(start)
            motion.phase{2,p}=-1;
        else
             motion.phase{2,p}=1;
        end
        % evtl noch 0 einfügen!
    end
    
    
%     motion_minV = abs(motion_min*ones(size(motion.r)));    
%     motion.r = motion.r + motion_minV;
%     
%     for p=1:NumOfPhases   
%     motion.phase{1,p} = motion.phase{1,p} + abs(motion_min*ones(1,2))   ;
%     end
    
    r1 = motion.r(1);
    r2 = motion.r(2);
    help = 0;
    r = motion.r(1);
    direction = 0;
    if(r1 < r2)
        for p = 1:NumOfPhases
            if (r > motion.phase{1,p}(1) && r < motion.phase{1,p}(2))
                motion.p(1) = p;
                help = 1;
                direction = 1;
            end
        end
    else
        for p = 1:NumOfPhases
            if (r < motion.phase{1,p}(1) && r > motion.phase{1,p}(2))
                motion.p(1) = p;
                help = 1;
                direction = -1;
            end
        end
    end
    if(help == 0)
        motion.p(1) = 1;
    end
    
    j = 2;
    p=  motion.p(1);
    while (j < length(motion.time))
        if( direction == -1)
        while(motion.r(j) > motion.phase{1,p}(2) && motion.r(j) > motion.r(j-1) )
            motion.p(j) = p;
            j = j+1;
        end
        p = p+1;
        end
        %HIER FEHLEN NOCH VIELE FÄLLE!!!
    end
       
end % if motion = CALYPSO;
    
    
    
    %%%%%%%%%%%%%%%END MOTION
    
    
    
    
end %beams


%TestTestTestTestTest
% for i=1:length(delivery)
%     for j = 1: length(delivery(i).time)
%         delivery(i).phase(j) = 1;
%     end
% end
%TestTestTestTestTest



%berechne Dosis für jede Phase  
disp('calc dose in each CT phase')
for p=1:NumOfPhases
    w=zeros(dij.totalNumOfBixels, 1); 
    for i=1:length(delivery)
        t =  find(delivery(i).phase == p  & ~isnan(delivery(i).w));
             
        w(delivery(i).j(t)) = delivery(i).w(t);
        Fluenzw{p} = w;
    end
 
    %Test    
    %resultGUI.phaseDose{p} = reshape(dij.physicalDose{1,1} * w, dij.dimensions);
    resultGUI.phaseDose{p} = reshape(dij.physicalDose{p} * w, dij.dimensions);
   
end

%just for testing: dose summation
% resultGUI.sumDose = resultGUI.phaseDose{1}(:);
% for p=2:NumOfPhases
%     resultGUI.sumDose = resultGUI.sumDose + resultGUI.phaseDose{p}(:);
% end
% 
% resultGUI.sumDose = reshape(resultGUI.sumDose, dij.dimensions);
    


