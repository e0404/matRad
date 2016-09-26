function [resultGUI, delivery, ct] = matRad_calc4dDose(ct, FileName)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad 4D dose calculation
% 
% call
%   
%
% input 
%       FileName = Name of PBP and Lmdout File
%   
%
% output
%   
%
% References
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~isdeployed % only if _not_ running as standalone
    
    % add path for optimization functions
    matRadRootDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(matRadRootDir,'4Ddose'))
    
    addpath(fullfile(matRadRootDir,'tools'))
    
    % get handle to Matlab command window
    mde         = com.mathworks.mde.desk.MLDesktop.getInstance;
    cw          = mde.getClient('Command Window');
    xCmdWndView = cw.getComponent(0).getViewport.getComponent(0);
    h_cw        = handle(xCmdWndView,'CallbackProperties');

    % set Key Pressed Callback of Matlab command window
    set(h_cw, 'KeyPressedCallback', @matRad_CWKeyPressedCallback);

end


%reads in PB XML Plan and result of dose delivery simulation and creates
%delievery struct
delivery = matRad_readLmdout(FileName);

%dose in each CT phase is calculated
[resultGUI, delivery] = matRad_calcPhaseDose(delivery);

%dose accumulation
[dAcc, ct] = matRad_doseAcc(ct, resultGUI,'DDM');  %acc Methods: 'EMT' 'DDM'
resultGUI.accDose = dAcc;

%visualisation
matRad_plotPhaseDose(resultGUI, 110); %TKUH005 slice 110 % T6H slice 50








