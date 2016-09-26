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


   

    %Motion assumption
    NumOfPhases = multScen.numOfCtScen;
    %Annahme 5s und linear  ÄNDERN
    phaseTime = 3/NumOfPhases;
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


%TestTestTestTestTest
for i=1:length(delivery)
    for j = 1: length(delivery(i).time)
        delivery(i).phase(j) = 1;
    end
end
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
    resultGUI.phaseDose{p} = reshape(dij.physicalDose{p,1} * w, dij.dimensions);
   
end

%just for testing: dose summation
resultGUI.sumDose = resultGUI.phaseDose{1}(:);
for p=2:NumOfPhases
    resultGUI.sumDose = resultGUI.sumDose + resultGUI.phaseDose{p}(:);
end

resultGUI.sumDose = reshape(resultGUI.sumDose, dij.dimensions);
    


