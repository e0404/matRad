% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrad_simulationAutoRunDaemon enables continuos evaluation of simulation
% 
% call
%   matrad_simulationAutoRunDaemon
%
% usage
%   place this file in the parent folder of your simulations and create a 
%   directory 'queue'. Place patients .mat in a subfolder as well as a 
%   setupStudy configuration file. This daemon will continue to look for
%   new simulation folders until the script is terminated by the user.
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

directoryOfQueue = 'queue';
marker = {'_RUNNING', '_FINISHED', '_ERROR'};
numWorker = 4;

currPath = pwd;
computeQueue = true;

while computeQueue
    listOfPatients = dir(directoryOfQueue);
    listOfPatients(1:2) = [];  % remove '.' and '..'
    
    % if no patient is in queue wait a bit longer
    if numel(listOfPatients) == 0
        fprintf('No patient found. Waiting ... \n');
        pause(300);
        continue
    end
    
    % reduce list to non-running and non-finished
    i = 1;
    while true
        if (~isempty(strfind(listOfPatients(i).name, marker{1})) ...
                || ~isempty(strfind(listOfPatients(i).name, marker{2})) ...
                || ~isempty(strfind(listOfPatients(i).name, marker{3})))
            listOfPatients(i) = [];
            i = 1; % jump to beginning of list
        else
            i = i + 1; % step forward
        end
        if i >= numel(listOfPatients) + 1
            break
        end
    end
    
    % if all patients are running or finished
    if numel(listOfPatients) == 0
        fprintf('All simulations finished or already running. Waiting ... \n');
        pause(300);
        continue
    end
    
    % pick a patient randomly
    currPatientIx = randi(numel(listOfPatients));
    currPatientPath = fullfile(listOfPatients(currPatientIx).folder, listOfPatients(currPatientIx).name);
    tempPatientPath = [currPatientPath marker{1}];
    movefile(currPatientPath, tempPatientPath);
    
    % get or start parpool
    p = gcp('nocreate');
    if numel(p) == 0
        try
            parpool(numWorker);
        catch
            warning('Could not start parpool.');
        end
    end
    
    cd(tempPatientPath)
    try
        setupStudy
        % after completion change directory back and rename
        cd(currPath);
        finishedPatientPath = [currPatientPath marker{2}];
        movefile(tempPatientPath, finishedPatientPath);
    catch
        warning('Error during computation of simulation.');
        cd(currPath);
        finishedPatientPath = [currPatientPath '_ERROR'];
        movefile(tempPatientPath, finishedPatientPath);
    end
    
    % safety wait to not overdo it
    pause(10)
end
