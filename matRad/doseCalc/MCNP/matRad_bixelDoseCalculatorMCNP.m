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

runDir = fullfile(this.workingDir, 'runfiles');
runFileList = dir(fullfile(runDir, 'MCNPrunfile_*bixel'));

switch this.externalCalculation
    case 'off'
        %% Run MCNP for every bixel in the run file directory
        if ~this.MCNPinstallationCheck
            matRad_cfg.dispError(['MCNP simulation requested but no MCNP installation found on this computer! ' ...
                'Set externalCalculation to ''write'' to only generate the run files.']);
        end

        oldDir = cd(runDir);
        restoreDir = onCleanup(@() cd(oldDir));

        parfor bixelCounter = 1:size(runFileList,1)
            matRad_cfg.dispInfo('MCNP calculation of dose distribution for bixel %d of %d...\n', bixelCounter, size(runFileList,1));
            tic;
            if ispc
                system(['mcnp6 I=', runFileList(bixelCounter).name, ...
                    ' OUTP=', runFileList(bixelCounter).name, 'o ', ...
                    ' RUNTPE=', runFileList(bixelCounter).name, 'r ', ...
                    ' MCTAL=', runFileList(bixelCounter).name, 'm ', ...
                    'MDATA= ', runFileList(bixelCounter).name, 'd']);
            else
                system(['mcnp6 I=', runFileList(bixelCounter).name, ...
                    ' OUTP=', runFileList(bixelCounter).name, 'o ', ...
                    ' RUNTPE=', runFileList(bixelCounter).name, 'r ', ...
                    ' MCTAL=', runFileList(bixelCounter).name, 'm ', ...
                    ' MESHTAL=', runFileList(bixelCounter).name, 'meshtal ']);
            end
            delete(strcat(runFileList(bixelCounter).name, 'o'));
            delete(strcat(runFileList(bixelCounter).name, 'r'));
            delete(strcat(runFileList(bixelCounter).name, 'd'));
            matRad_cfg.dispInfo('Calculation for bixel %d took %g seconds.\n', bixelCounter, toc);
        end

    case 'write'
        %% Only write a script to run all bixels externally
        cores = feature('numcores');    % Attention: should be adopted to allow portability to other pc/cluster
        if ispc
            scriptName = fullfile(runDir, 'runAll.cmd');
        else
            scriptName = fullfile(runDir, 'runAll.sh');
        end
        commandMCNP = 'mpiexec -np %d mcnp6.mpi n=MCNPrunfile_%dbixel\n';
        fileID_runAll = fopen(scriptName, 'w');
        for i = 1:size(runFileList,1)
            fprintf(fileID_runAll, commandMCNP, cores, i);
        end
        fclose(fileID_runAll);
        matRad_cfg.dispInfo('MCNP simulation skipped for external calculation.\nRun files and %s have been written to: "%s"\n', scriptName, strrep(runDir,'\','\\'));

    otherwise
        matRad_cfg.dispError('Unknown externalCalculation mode ''%s''! Use ''off'', ''write'' or a folder with existing results.', this.externalCalculation);
end
end
