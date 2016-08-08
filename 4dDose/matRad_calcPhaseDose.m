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

% get data from workspace
multScen       = evalin('base','multScen');
resultGUI      = evalin('base','resultGUI');
dij            = evalin('base','dij');

for i=1:length(delivery)
    delivery(i).phase = zeros(length(delivery(i).es), 1); 
    delivery(i).w = zeros(length(delivery(i).es), 1); 
    delivery(i).index_j = zeros(length(delivery(i).es), 1); 
    
    %für richtige index der Bixel in Fluenzvektor
    add = 0;
    if(i>1)        
        for k=2:i
            add = length(delivery(k-1).w);       
        end
    end
    
    [j, index_j] = sort(delivery(i).j);
    delivery(i).w = resultGUI.finalw(index_j + add);
    delivery(i).index_j = index_j + add;

  

    %Motion assumption
    NumOfPhases = multScen.numOfCtScen;
    %Annahme 5s und linear  ÄNDERN
    phaseTime = 5/NumOfPhases;
    j=1;
    p = 1;
    Time = 0;
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
end %beams

%berechne Dosis für jede Phase  
disp('calc dose in each CT phase')
for p=1:NumOfPhases
    w=zeros(dij.totalNumOfBixels, 1); 
    for i=1:length(delivery)
        t =  find(delivery(i).phase == p);
        w(delivery(i).index_j(t)) = delivery(i).w(t);
        Fluenzw{p} = w;
    end
 
    %Test    
    %resultGUI.phaseDose{p} = reshape(dij.physicalDose{1,1} * w, dij.dimensions);
    resultGUI.phaseDose{p} = reshape(dij.physicalDose{p,1} * w, dij.dimensions);
   
end

%just for testing: dose summation
resultGUI.sumDose = resultGUI.phaseDose{1}(:);
for p=2:NumOfPhases
    resultGUI.sumDose = resultGUI.sumDose + resultGUI.phaseDose{p}(:);
end

resultGUI.sumDose = reshape(resultGUI.sumDose, dij.dimensions);
    


