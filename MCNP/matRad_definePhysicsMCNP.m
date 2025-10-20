function matRad_definePhysicsMCNP(fileID_C_rest, this, binIntervals, simPropMCNP)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define physics that is needed in MCNP simulation according to the used
% tally and the radiation type.
%
% Neutron dose engine A: Monte Carlo - MCNP6
% 
% call
%   matRad_definePhysicsMCNP(fileID_C_rest)
%
% input
%   fileID_C_rest:  File ID used for block C
%
% output
%
%   none
%
% References
%   [1] PELOWITZ, D. B., et al. MCNP6 User’s Manual. LACP-00634, May, 2013.
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 12/2018
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set Mode and Phys cards for particle transport
% Set up TMESH tally for total dose calculation
        % if isfield(this.propOpt, 'bioOptimization')
        %     fprintf(fileID_C_rest, 'MODE N P E H D T S A #\n');
        %         fprintf(fileID_C_rest, 'PHYS:N 100 0 0 J J J 4 -1 1 J J 0 0\n');
        %         fprintf(fileID_C_rest, 'PHYS:P 100 0 0 1 0 J 0\n');
        %         fprintf(fileID_C_rest, 'PHYS:E 100 0 0 0 0 1 1 1 1 0 0 J J 0.917 0.001 0\n');
        %         fprintf(fileID_C_rest, 'PHYS:H 100 0 -1 J 0 J 1 3J 0 0 0 0.917 J 0 -1\n');
        %         fprintf(fileID_C_rest, 'PHYS:D 16J -1\n');  % change default settings s. th. delta-rays production is on
        %         fprintf(fileID_C_rest, 'PHYS:T 16J -1\n');
        %         fprintf(fileID_C_rest, 'PHYS:S 16J -1\n');
        %         fprintf(fileID_C_rest, 'PHYS:A 16J -1\n');
        %         fprintf(fileID_C_rest, 'PHYS:# 16J -1\n');
        %         fprintf(fileID_C_rest, 'CUT:N J 0.00000000001\n');
        %         fprintf(fileID_C_rest, 'CUT:P J 0.003\n');
        %         fprintf(fileID_C_rest, 'CUT:E J 0.00001\n');
        %         fprintf(fileID_C_rest, 'CUT:H J 0.0009\n');
        %         fprintf(fileID_C_rest, 'CUT:D J 0.0009\n');
        %         fprintf(fileID_C_rest, 'CUT:T J 0.0009\n');
        %         fprintf(fileID_C_rest, 'CUT:S J 0.0009\n');
        %         fprintf(fileID_C_rest, 'CUT:A J 0.0009\n');
        %         fprintf(fileID_C_rest, 'CUT:# J 0.0009\n');           
        % else
    fprintf(fileID_C_rest, 'MODE N P E H D T S A #\n');
    fprintf(fileID_C_rest, 'PHYS:N 100 0 0 J J J 4 -1 1 J J 0 0\n');
    fprintf(fileID_C_rest, 'PHYS:P 100 0 0 1 0 J 0\n');
    fprintf(fileID_C_rest, 'PHYS:E 100 0 0 0 0 1 1 1 1 0 0 J J 0.917 0.001 0\n');
    fprintf(fileID_C_rest, 'PHYS:H 100 0 -1 J 0 J 1 3J 0 0 0 0.917 J 0 -1\n');
    fprintf(fileID_C_rest, 'PHYS:D 16J -1\n');  % change default settings s. th. delta-rays production is on
    fprintf(fileID_C_rest, 'PHYS:T 16J -1\n');
    fprintf(fileID_C_rest, 'PHYS:S 16J -1\n');
    fprintf(fileID_C_rest, 'PHYS:A 16J -1\n');
    fprintf(fileID_C_rest, 'PHYS:# 16J -1\n');
    fprintf(fileID_C_rest, 'CUT:N J 0.00000000001\n');
    fprintf(fileID_C_rest, 'CUT:P J 0.003\n');
    fprintf(fileID_C_rest, 'CUT:E J 0.05\n');
    fprintf(fileID_C_rest, 'CUT:H J 1\n');
    fprintf(fileID_C_rest, 'CUT:D J 1\n');
    fprintf(fileID_C_rest, 'CUT:T J 1\n');
    fprintf(fileID_C_rest, 'CUT:S J 1\n');
    fprintf(fileID_C_rest, 'CUT:A J 1\n');
    fprintf(fileID_C_rest, 'CUT:# J 1\n');
%     case 'KERMA_F4'
%         switch pln.radiationMode
%             case 'neutrons'
%                 fprintf(fileID_C_rest, 'MODE N\n');
%                 fprintf(fileID_C_rest, 'PHYS:N 100 0 0 J J J 0 -1 0 J J 1 0\n');
%                 fprintf(fileID_C_rest, 'CUT:N J 0.000000015\n');              
%             case 'photons'
%                 fprintf(fileID_C_rest, 'MODE P\n');
%                 fprintf(fileID_C_rest, 'PHYS:P 100 0 0 0 0 J 0\n');
%                 fprintf(fileID_C_rest, 'CUT:P J 0.005\n');
%             case'neutrons_PLUS_photons'
%                 fprintf(fileID_C_rest, 'MODE N P\n');
%                 fprintf(fileID_C_rest, 'PHYS:N 100 0 0 J J J 0 -1 1 J J 0 0\n');
%                 fprintf(fileID_C_rest, 'PHYS:P 100 0 0 0 0 J 0\n');
%                 fprintf(fileID_C_rest, 'CUT:N J 0.000000015\n');   
%                 fprintf(fileID_C_rest, 'CUT:P J 0.005\n');                
%         end
% end

% Set random number generator properties
fprintf(fileID_C_rest, 'RAND GEN=2 SEED=43 STRIDE=10000000001\n');   % Gen 2: period= 9.2e18 
                                                                 % Seed makes sense
                                                                 % Stride: Number of random numbers between source particles
switch simPropMCNP.geometryOption
    case 'Lattice'                                                                                               
        % Define temperature to be at body temperature
        fprintf(fileID_C_rest, 'TMP 2.53e-08\n');   % air voxels are set to room temp
        for iCellTemp=2:size(binIntervals,2)
            fprintf(fileID_C_rest, '        0.0000000266911575\n');   % tissue voxels are set to body temp = 309.75 K in MeV
        end
        for iCellImportance=size(binIntervals,2):size(binIntervals,2)+2
            fprintf(fileID_C_rest, '        2.53e-08\n');   % set temperature for helper cells
        end
            fprintf(fileID_C_rest, '        2.53e-08\n');   % set temperature for cell outside sim volume
        
        
        % Define importance for each assigned material individually
        fprintf(fileID_C_rest, 'IMP:N,P,E,H,D,T,S,A,# \n');   % same importance for all particles
        for iCellImportance=1:size(binIntervals,2)
            fprintf(fileID_C_rest, ['        ',num2str(binIntervals(iCellImportance).importance), '\n']);   % set importance for all particle for each material type individually
        end
        for iCellImportance=size(binIntervals,2):size(binIntervals,2)+2
            fprintf(fileID_C_rest, ['        ','1', '\n']);   % set importance in helper cells to 1
        end
        fprintf(fileID_C_rest, ['        ','0', '\n']);   % set importance to zero outside sim volume

        
    otherwise
        error('Geometry option not supported.')
end

fprintf(fileID_C_rest, ['NPS' ' ' num2str(this.config.Num_Primaries) '\n']);