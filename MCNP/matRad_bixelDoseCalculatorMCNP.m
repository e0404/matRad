function matRad_bixelDoseCalculatorMCNP(this)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad neutron dose calculation for each bixel individually
%
% Neutron dose engine A: Monte Carlo - MCNP6
%
% call
%   dij = matRad_calcPhotonDose(pathToRunfiles, stf, ct, pln, cst, binIntervals)
%
% input
%   pathToRunfiles: indicate path to MCNP runfiles here
%   stf, ct, pln, binIntervals
%
% output
%   dij:            matRad dij struct
%
% References
%   [1] PELOWITZ, D. B., et al. MCNP6 Users Manual. LACP-00634, May, 2013.
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 11/2018
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matRad_cfg = MatRad_Config.instance();

if this.MCNPinstallationCheck && ~this.externalCalculation
    %% Go to runfiles and get list of runfiles within directory
    cd(strcat(matRad_cfg.matRadRoot, filesep, 'MCNP', filesep, 'runfiles_tmp'));
    runFileList = dir('MCNPrunfile_*bixel');

    wb = waitbar(0, ['Calculating dose for bixel: ', num2str(1)], 'Name', 'Dose Calculation with MCNP');

    %% Run calculation for every bixel
    parfor bixelCounter=1:size(runFileList,1)

        disp('*****')
        disp(['MCNP calculation of dose distribution for bixel ', num2str(bixelCounter), '...'])
        disp('*****')

        tic;

        %     waitbar(bixelCounter/size(runFileList,1), wb, ['Calculating dose for bixel: ', num2str(bixelCounter)], 'Name', 'Dose Calculation with MCNP');
        if ispc 
            % system(['mpiexec -np ',num2str(numberCores4U), ' mcnp6.mpi I=', runFileList(bixelCounter).name, ...
            system(['mcnp6 I=', runFileList(bixelCounter).name, ...
                 ' OUTP=', runFileList(bixelCounter).name, 'o ', ...
                 ' RUNTPE=', runFileList(bixelCounter).name, 'r ', ...
                 ' MCTAL=', runFileList(bixelCounter).name, 'm ', ...
                 'MDATA= ', runFileList(bixelCounter).name, 'd']);
            % Clean up
            delete(strcat(runFileList(bixelCounter).name, 'o'))
            delete(strcat(runFileList(bixelCounter).name, 'r'))
            delete(strcat(runFileList(bixelCounter).name, 'd'))

        else   
            system(['mcnp6 I=', runFileList(bixelCounter).name, ...
                ' OUTP=', runFileList(bixelCounter).name, 'o ', ...
                ' RUNTPE=', runFileList(bixelCounter).name, 'r ', ...
                ' MCTAL=', runFileList(bixelCounter).name, 'm ', ...
                ' MESHTAL=', runFileList(bixelCounter).name, 'meshtal '])
            % Clean up
            delete(strcat(runFileList(bixelCounter).name, 'o'))
            delete(strcat(runFileList(bixelCounter).name, 'r'));
            delete(strcat(runFileList(bixelCounter).name, 'd'));
        end

        calculationTime = toc;
        disp('*****')
        disp(['Calculation for bixel ', num2str(bixelCounter), ' took ', num2str(calculationTime), ' seconds.'])
        disp('*****')

    end

    close(wb)
elseif this.externalCalculation
    if ispc % Write script to run MCNP simulation
        cores = feature('numcores');    % Attention: should be adopted to allow portability to other pc/cluster
        cd(strcat(matRad_cfg.matRadRoot, filesep, 'MCNP', filesep, 'runfiles_tmp'));
        runFileList = dir('MCNPrunfile_*bixel');
        fileID_runAll = fopen('runAll.cmd', 'w');
        for i=1:size(runFileList,1); fprintf(fileID_runAll, strcat('mpiexec -np ',in2str(cores),' mcnp6.mpi n=MCNPrunfile_bixel', int2str(i), '\n')); end
        fclose(fileID_runAll);
    end

    matRad_cfg.dispInfo('Please use question dialog to continue after finishing external calculations.\n')
    matRad_cfg.dispInfo('*****\n')
    % External calculation
    answer = questdlg('Did external MCNP simulations finish?', ...
        'External Calculation', ...
        'Yes', 'No', 'No');

    while ~strcmp(answer, 'Yes')
        matRad_cfg.dispInfo('matRad will crash if you continue without finishing external calculation.\n')
        answer = questdlg('Did external MCNP simulations finish?', ...
            'External Calculation', ...
            'Yes', 'No', 'No');
    end
elseif ~this.externalCalculation && this.MCNPinstallationCheck
    matRad_cfg.dispWarning('MCNP simulation requested but no MCNP installation found on your computer!\n')
end
end