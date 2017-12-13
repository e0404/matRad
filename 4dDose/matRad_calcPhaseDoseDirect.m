function [resultGUI, delivery] = matRad_calcPhaseDoseDirect(ct, stf, pln, cst, resultGUI, dij, delivery)
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


    
    
        NumOfPhases = delivery.NumOfPhases; %ct.numOfCtScen;
        
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
    
  
       
end


numOfCtScen = ct.numOfCtScen;

help = 0;
if(numOfCtScen ~= NumOfPhases)
    help = 1;
end

%berechne Dosis für jede Phase  
disp('calc dose in each CT phase')
pln.multScen.numOfCtScen = 1; 
ctPhase.numOfCtScen = 1;
ctPhase.resolution = ct.resolution;
ctPhase.cubeDim = ct.cubeDim;
for p=1:NumOfPhases
      w=zeros(dij.totalNumOfBixels, 1); 
      for i=1:length(delivery)
         t =  find(delivery(i).phase(:,p) ~= 0  & ~isnan(delivery(i).w));
         w(delivery(i).j(t)) = delivery(i).w(t) .* delivery(i).phase(t,p);       
      end
 
      ctPhase.cube{1} = ct.cube{1,p+help};
      disp(['Dose calculation in CT phase ', num2str(p)]);
       resultHelp = matRad_calcDoseDirect(ctPhase,stf,pln,cst,w);
       resultGUI.phaseDose{p} = resultHelp.physicalDose;
end
  
pln.multScen.numOfCtScen = numOfCtScen;
ctPhase.numOfCtScen = numOfCtScen;

