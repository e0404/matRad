function [resultGUI, delivery] = matRad_calcPhaseDoseMatrix(resultGUI, dij, delivery, MOTION, phaseTimeDev)
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
% Silke Ulrich May 2017
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MOTION =  'linear';  %'sampled_periode' ; %'linear'; %'sampled_phase' 
    

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
        NumOfPhases = delivery(1).NumOfPhases; %size(dij.physicalDose,1);
        offset = delivery(1).offset; %0;
        phaseTime = delivery(1).motionperiod/NumOfPhases; %5/NumOfPhases;
         
        delivery(i).phase = zeros(length(delivery(i).time),NumOfPhases);
    
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
        if (delivery(i).time(t)/1000000 < Time)  
            while(t < length(delivery(i).time) && delivery(i).time(t)/1000000 < Time)  
            delivery(i).phase(t, p) = 1;
            t= t+1;
            end
            Time = Time + phaseTime;
            if(t < length(delivery(i).time) && delivery(i).time(t)/1000000 < Time)
                delivery(i).phase(t, p) = (Time - phaseTime - delivery(i).time(t-1)/1000000)/(delivery(i).time(t)/1000000 - delivery(i).time(t-1)/1000000);
                p = p+1;
                if(p > NumOfPhases)
                p=1;
                end 
                delivery(i).phase(t, p) = (delivery(i).time(t)/1000000 - (Time - phaseTime))/(delivery(i).time(t)/1000000 - delivery(i).time(t-1)/1000000);
                t = t+1;
            else
                p = p+1;
                if(p > NumOfPhases)
                    p=1;
                end 
            end
        else
            p = p+1;
            if(p > NumOfPhases)
                    p=1;
                end 
            Time = Time + phaseTime;
        end
    end
    
    elseif(strcmp(MOTION,'sampled_periode')) 
        NumOfPhases = delivery(1).NumOfPhases;
     % NumOfPhases = size(dij.physicalDose,1);
        offset = delivery(1).offset; %0;
        
        phaseTimeMean = delivery(1).motionperiod/NumOfPhases; 
        %phaseTimeDev = 0.05; %???
    
        phaseTime = (phaseTimeDev .*randn(1) + phaseTimeMean);
         
        delivery(i).phase = zeros(length(delivery(i).time),NumOfPhases);
        t=1;
        p = ceil((offset+0.001)/phaseTime);  
        while p > NumOfPhases
            p = p-NumOfPhases;
        end
    
        Time = mod(offset, phaseTime); 
        if(Time == 0)
            Time = Time + phaseTime;
        end
        while(t < length(delivery(i).time))       
            if (delivery(i).time(t)/1000000 < Time)  
                while(t < length(delivery(i).time) && delivery(i).time(t)/1000000 < Time)  
                    delivery(i).phase(t, p) = 1;
                    t= t+1;
                end
                %Time = Time + phaseTime;
                if(t < length(delivery(i).time) && delivery(i).time(t)/1000000 < Time+phaseTime)
                    delivery(i).phase(t, p) = (Time - delivery(i).time(t-1)/1000000)/(delivery(i).time(t)/1000000 - delivery(i).time(t-1)/1000000);
                    p = p+1;
                    if(p > NumOfPhases)
                        p=1;
                        phaseTime = phaseTimeDev .*randn(1) + phaseTimeMean;
                    end 
                    delivery(i).phase(t, p) = (delivery(i).time(t)/1000000 - (Time))/(delivery(i).time(t)/1000000 - delivery(i).time(t-1)/1000000);
                    t = t+1;
                else
                    p = p+1;
                    if(p > NumOfPhases)
                        p=1;
                        phaseTime = phaseTimeDev .*randn(1) + phaseTimeMean;
                    end 
                end
                
            else
                p = p+1;
                if(p > NumOfPhases)
                    p=1;
                    phaseTime = phaseTimeDev .*randn(1) + phaseTimeMean;
                end 
            end
            Time = Time + phaseTime;
        end        
   
    
    elseif(strcmp(MOTION,'sampled_phase'))
        NumOfPhases = delivery(1).NumOfPhases;
      %NumOfPhases = size(dij.physicalDose,1);
        offset = delivery(1).offset; %0;
        phaseTimeMean = delivery(1).motionperiod/NumOfPhases; 
        %phaseTimeDev = 0.05; %?????
    
        %phaseTime =  (phaseTimeDev .*randn(1) + phaseTimeMean)
        delivery(i).phase = zeros(length(delivery(i).time),NumOfPhases);
    
        t=1;
        p = ceil((offset+0.001)/(phaseTimeDev .*randn(1) + phaseTimeMean));  
        while p > NumOfPhases
            p = p-NumOfPhases;
        end
    
        Time = mod(offset, (phaseTimeDev .*randn(1) + phaseTimeMean)); 
        if(Time == 0)
            Time = Time + (phaseTimeDev .*randn(1) + phaseTimeMean);
        end
        while(t < length(delivery(i).time))       
            if (delivery(i).time(t)/1000000 < Time)  
                while(t < length(delivery(i).time) && delivery(i).time(t)/1000000 < Time)  
                delivery(i).phase(t, p) = 1;
                t= t+1;
                end
                %Time = Time + (phaseTimeDev .*randn(1) + phaseTimeMean);
                if(t < length(delivery(i).time) && delivery(i).time(t)/1000000 < Time + (phaseTimeDev .*randn(1) + phaseTimeMean))
                    delivery(i).phase(t, p) = (Time - delivery(i).time(t-1)/1000000)/(delivery(i).time(t)/1000000 - delivery(i).time(t-1)/1000000);
                    p = p+1;
                    if(p > NumOfPhases)
                        p=1;
                    end 
                    delivery(i).phase(t, p) = (delivery(i).time(t)/1000000 - Time )/(delivery(i).time(t)/1000000 - delivery(i).time(t-1)/1000000);
                    t = t+1;
                else
                    p = p+1;
                    if(p > NumOfPhases)
                        p=1;
                    end 
                end
            else
                p = p+1;
                if(p > NumOfPhases)
                    p=1;
                end                 
            end
            Time = Time + (phaseTimeDev .*randn(1) + phaseTimeMean);
        end  
    end   
end

    
%berechne Dosis für jede Phase  
%disp('calc dose in each CT phase')
help = 0;
if(size(dij.physicalDose, 1) ~= NumOfPhases)
    help = 1;
end
for p=1:NumOfPhases
      w=zeros(dij.totalNumOfBixels, 1); 
      for i=1:length(delivery)
         t =  find(delivery(i).phase(:,p) ~= 0  & ~isnan(delivery(i).w));
         w(delivery(i).j(t)) = delivery(i).w(t) .* delivery(i).phase(t,p);       
      end
 
       resultGUI.phaseDose{p} = reshape(dij.physicalDose{p+help} * w , dij.dimensions);
          
      
    if isequal(resultGUI.bioParam.type,'MCN_RBExD')
    
        resultGUI.phaseAlphaDose{p} = reshape(dij.mAlphaDose{p+help} * w, dij.dimensions);
        resultGUI.phaseSqrtBetaDose{p} = reshape(dij.mSqrtBetaDose{p+help} * w, dij.dimensions);
    
    end
end
  


