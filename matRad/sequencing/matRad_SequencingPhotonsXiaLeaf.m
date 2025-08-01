classdef  matRad_SequencingPhotonsXiaLeaf < matRad_SequencingPhotonsAbstract

% multileaf collimator leaf sequence algorithm 
% for intensity modulated beams with multiple static segments according to 
% Xia et al. (1998) Medical Physics
% 
%   [1] http://online.medphys.org/resource/1/mphya6/v25/i8/p1424_s1
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (Constant)
        name = 'Photons Xia Leaf Sequenceer';
        shortName = 'xia';
        possibleRadiationModes = {'photons'};
    end 
    properties
        mode = 'rl' % sliding window (sw) or reducing level (rl)
    end

    methods

        function sequence = sequence(this,w,stf)
            matRad_cfg = MatRad_Config.instance();
            numOfBeams = numel(stf);
            offset = 0;

            for i = 1:numOfBeams

                [D_0,D_k, shapes,calFac,indInMx] = this.initBeam(stf(i),w(1+offset:stf(i).numOfRays+offset));

                
                % Save the maximun intensity (Equation 5)
                L_k = max(D_k(:));
                
                % Save the maximun initial intensity matrix value in L_0.
                L_0 = L_k;
                
                % Set k=0, this variable is used for residuals intensity matrices D_k.
                k = 0;
                
                % start sequencer
                while L_k > 0
                    
                    k = k + 1;
                    
                    
                    %Rounded off integer. Equation 7.
                    m = floor(log2(L_k));
                    
                    % Convert m=1 if is less than 1. This happens when L_k belong to ]0,2[
                    if m < 1
                        m = 1;
                    end
                    
                    %Calculate the delivery intensity unit. Equation 6.
                    d_k = floor(2^(m-1));
                    
                    % Opening matrix.
                    openingMx = D_k >= d_k;
                    
                    dimOfFluenceMxZ = size(shapes,1);

                    switch this.mode
                        case 'sw' % sliding window technique!
                            for j = 1:dimOfFluenceMxZ
                                openIx = find(openingMx(j,:) == 1,1,'first');
                                if ~isempty(openIx)
                                    closeIx = find(openingMx(j,openIx+1:end) == 0,1,'first');
                                    if ~isempty(closeIx)
                                        openingMx(j,openIx+closeIx:end) = 0;
                                    end
                                end
                                
                            end
                       case  'rl' % reducing levels technique!
                            for j = 1:dimOfFluenceMxZ
                                [maxVal,maxIx] = max(openingMx(j,:) .* D_k(j,:));
                                if maxVal > 0
                                    closeIx = maxIx + find(openingMx(j,maxIx+1:end) == 0,1,'first');
                                    if ~isempty(closeIx)
                                        openingMx(j,closeIx:end) = 0;
                                    end
                                    openIx = find(openingMx(j,1:maxIx-1) == 0,1,'last');
                                    if ~isempty(openIx)
                                        openingMx(j,1:openIx) = 0;
                                    end
                                end
                                                
                            end
                       otherwise
                            matRad_cfg.dispError('unknown sequencing mode')
                       end
                                            
                    shape_k       = openingMx * d_k;
                                       
                    shapes(:,:,k) = shape_k;
                    shapesWeight(k) = d_k;
                    D_k = D_k - shape_k;
                    
                    L_k = max(D_k(:)); % eq 5
                    
                end
                
                sequence.beam(i).numOfShapes  = k;
                sequence.beam(i).shapes       = shapes(:,:,1:k);
                sequence.beam(i).shapesWeight = shapesWeight(1:k)/this.sequencingLevel*calFac;
                sequence.beam(i).bixelIx      = 1+offset:stf(i).numOfRays+offset;
                sequence.beam(i).fluence      = D_0;
                
                sequence.w(1+offset:stf(i).numOfRays+offset,1) = D_0(indInMx)/this.sequencingLevel*calFac;
            
                offset = offset + stf(i).numOfRays;
            
            end
            if this.visMode
               this.plotSegments(sequence)
            end
        end
    end
    methods  (Static)
        function [available,msg] = isAvailable(pln,machine)
            % see superclass for information            
                   
            if nargin < 2
                machine = matRad_loadMachine(pln);
            end

            % Check superclass availability
            [available,msg] = matRad_SequencingPhotonsAbstract.isAvailable(pln,machine);

            if ~available
                return;
            else
                available = false;
                msg = [];
            end
    
            %checkBasic
            try
                checkBasic = isfield(machine,'meta') && isfield(machine,'data');
    
                %check modality
                checkModality = any(strcmp(matRad_SequencingPhotonsXiaLeaf.possibleRadiationModes, machine.meta.radiationMode)) && any(strcmp(matRad_SequencingPhotonsXiaLeaf.possibleRadiationModes, pln.radiationMode));
                
                %Sanity check compatibility
                if checkModality
                    checkModality = strcmp(machine.meta.radiationMode,pln.radiationMode);
                end
    
                preCheck = checkBasic && checkModality;
    
                if ~preCheck
                    return;
                end
            catch
                msg = 'Your machine file is invalid and does not contain the basic field (meta/data/radiationMode)!';
                return;
            end

            available = preCheck;
        end
    end
end