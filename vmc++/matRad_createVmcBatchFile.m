function matRad_createVmcBatchFile(parallelSimulations,filepath,verboseLevel)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad batchfile creation
%
% call
%   matRad_createVmcBatchFile(parallelSimulations,filepath,verboseLevel)
%
% input
%   parallelSimulations: no of parallel simulations
%   filepath:            path where batchfile is created (this has to be the 
%                        path of the vmc++ folder)
%   verboseLevel:        optional. number specifying the amount of output 
%                        printed to the command prompt
%
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set verbose level
if nargin < 3
    verboseString = '/B'; % to not open terminals by default
else
    if verboseLevel > 1 % open terminals is verboseLevel > 1
        verboseString = '';
    else
        verboseString = '/B';
    end        
end

if ispc % parallelization only possible on windows systems

    parallelProcesses = cell(1,parallelSimulations);
    for i = 1:parallelSimulations
        parallelProcesses{1,i} = ['start "" 9>"%lock%',num2str(i),'" ' verboseString ' .\bin\vmc_Windows.exe MCpencilbeam_temp_',num2str(i)];
    end

    batchFile = {...
        ['@echo off'],...
        ['setlocal'],...
        ['set "lock=%temp%\wait%random%.lock"'],...
        [''],...
        parallelProcesses{:},...
        [''],...
        [':Wait for all processes to finish (wait until lock files are no longer locked)'],...
        ['1>nul 2>nul ping /n 2 ::1'],...
        ['for %%N in (',strjoin(arrayfun(@(x) num2str(x),(1:parallelSimulations),'UniformOutput',false),' '),') do ('],...
        ['  (call ) 9>"%lock%%%N" || goto :Wait'],...
        [') 2>nul'],...
        [''],...
        ['del "%lock%*"'],...
        ['']...
        %,['echo Done - ready to continue processing']
        };
    
elseif isunix
    
    batchFile = {'./bin/vmc_Linux.exe MCpencilbeam_temp_1'};

end

% write batch file
fid = fopen(filepath,'wt');
for i = 1 : length(batchFile)
  fprintf(fid,'%s\n',batchFile{i});
end
fclose(fid);

if isunix
   system(['chmod a+x ' filepath]);
end

end
