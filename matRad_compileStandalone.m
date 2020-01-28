%Create Standalone with and without runtime

if ismac
    suffix = 'Mac';
    % Mac platform
elseif isunix
    suffix = 'Linux';
    % Linux platform
elseif ispc
    suffix = 'Windows';
    % Windows platform
else
    error('Platform not supported')
end

applicationCompiler('-package',['matRad' suffix]);

disp('Done');


%% Files to be compiled separated by space
%fileNames="matRadGUI.fig matRadGUI.m";





%% Final command to run
%mcc -m matRadGUI.fig matRadGUI.m

% if isunix && ~ismac
%     fname 'linux64'