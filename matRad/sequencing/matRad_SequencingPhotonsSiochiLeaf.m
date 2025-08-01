classdef  matRad_SequencingPhotonsSiochiLeaf < matRad_SequencingPhotonsAbstract

    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties (Constant)
        name = 'Photons Siochi Leaf Sequenceer';
        shortName = 'siochi';
        possibleRadiationModes = {'photons'};

    end 

    methods


        function sequence = sequence(this,w,stf)

            
            offset = 0;
                        
            for i = 1:numel(stf)
                
                [D_0,D_k, shapes,calFac,indInMx] = this.initBeam(stf(i),w(1+offset:stf(i).numOfRays+offset));
                
                shapesWeight = zeros(10000,1);
                k = 0;
                              
                %Decompose the port, do rod pushing
                [tops, bases] = this.decomposePort(D_k);
                %Form segments
                [shapes,shapesWeight,k]=this.convertToSegments(shapes,shapesWeight,k,tops,bases);
                    
                sequence.beam(i).numOfShapes  = k;
                sequence.beam(i).shapes       = shapes(:,:,1:k);
                sequence.beam(i).shapesWeight = shapesWeight(1:k)/this.sequencingLevel*calFac;
                sequence.beam(i).bixelIx      = 1+offset:stf(i).numOfRays+offset;
                sequence.beam(i).fluence      = D_0;
                sequence.beam(i).sum          = zeros(size(D_0));
                
                for j = 1:k
                    sequence.beam(i).sum = sequence.beam(i).sum+sequence.beam(i).shapes(:,:,j)*sequence.beam(i).shapesWeight(j);
                end
                sequence.w(1+offset:stf(i).numOfRays+offset,1) = sequence.beam(i).sum(indInMx);
                
                offset = offset + stf(i).numOfRays;
                
            end

            if this.visMode
                this.plotSegments(sequence)
            end
        end


        function [tops, bases] = decomposePort(~,map)
            %Returns tops and bases of a fluence matrix "map" for Siochi leaf
            %sequencing algorithm (rod pushing part).  Accounts for collisions and
            %tongue and groove (Tng) effects.
            
            [dimZ,dimX] = size(map);
            map_nonZero = (map~=0);

            [D_k_Z, D_k_X] = ind2sub([dimZ,dimX], find(map_nonZero));
            minZ = min(D_k_Z);
            maxZ = max(D_k_Z);
            minX = min(D_k_X);
            maxX = max(D_k_X);
                
            tops = zeros(dimZ, dimX);
            bases = zeros(dimZ, dimX);
            
            for i = minX:maxX
                maxTop = -1;
                TnG = 1;
                for j = minZ:maxZ
                    if i == minX
                        bases(j,i) = 1;
                        tops(j,i) = bases(j,i)+map(j,i)-1;
                    else %assign trial base positions
                        if map(j,i) >= map(j,i-1) %current rod >= previous, match the bases
                            bases(j,i) = bases(j,i-1);
                            tops(j,i) = bases(j,i)+map(j,i)-1;
                        else %current rod <previous
                            if map(j,i) == 0 %rod length=0, put in in next slab after top of previous
                                bases(j,i) = tops(j,i-1)+1;
                                tops(j,i) = bases(j,i)-1;
                            else %rod length~=0, match tops
                                tops(j,i) = tops(j,i-1);
                                bases(j,i) = tops(j,i)-map(j,i)+1;
                            end
                        end
                    end
                    %determine which rod has the highest top in column
                    if tops(j,i) > maxTop
                        maxTop = tops(j,i);
                        maxRow = j;
                    end
                end
                
                %Correct for collision and tongue and groove error
                while(TnG)
                    %go from maxRow down checking for TnG.  This occurs when a shorter
                    %rod is "peeking over" a longer one in the direction transverse to
                    %the leaf motion.  To fix this, match either the tops or bases of
                    %the rods.
                    for j = (maxRow-1):-1:minZ
                        if map(j,i) < map(j+1,i)
                            if tops(j,i) > tops(j+1,i)
                                tops(j+1,i) = tops(j,i);
                                bases(j+1,i) = tops(j+1,i)-map(j+1,i)+1;
                            elseif bases(j,i) < bases(j+1,i)
                                bases(j,i) = bases(j+1,i);
                                tops(j,i) = bases(j,i)+map(j,i)-1;
                            end
                        else
                            if tops(j,i) < tops(j+1,i)
                                tops(j,i) = tops(j+1,i);
                                bases(j,i) = tops(j,i)-map(j,i)+1;
                            elseif bases(j,i) > bases(j+1,i)
                                bases(j+1,i) = bases(j,i);
                                tops(j+1,i) = bases(j+1,i)+map(j+1,i)-1;
                            end
                        end
                    end
                    %go from maxRow up checking for TnG
                    for j = (maxRow+1):maxZ
                        if map(j,i) < map(j-1,i)
                            if tops(j,i) > tops(j-1,i)
                                tops(j-1,i) = tops(j,i);
                                bases(j-1,i) = tops(j-1,i)-map(j-1,i)+1;
                            elseif bases(j,i) < bases(j-1,i)
                                bases(j,i) = bases(j-1,i);
                                tops(j,i) = bases(j,i)+map(j,i)-1;
                            end
                        else
                            if tops(j,i) < tops(j-1,i)
                                tops(j,i) = tops(j-1,i);
                                bases(j,i) = tops(j,i)-map(j,i)+1;
                            elseif bases(j,i) > bases(j-1,i)
                                bases(j-1,i) = bases(j,i);
                                tops(j-1,i) = bases(j-1,i)+map(j-1,i)-1;
                            end
                        end
                    end
                    %now check if all TnG conditions have been removed
                    TnG = 0;
                    for j = (minZ+1):maxZ
                        if map(j,i) < map(j-1,i);
                            if tops(j,i) > tops(j-1,i)
                                TnG = 1;
                            elseif bases(j,i) < bases(j-1,i)
                                TnG = 1;
                            end
                        else
                            if tops(j,i) < tops(j-1,i)
                                TnG = 1;
                            elseif bases(j,i) > bases(j-1,i)
                                TnG = 1;
                            end
                        end
                    end
                end
            end
        end

        function [shapes,shapesWeight,k] = convertToSegments(this, shapes,shapesWeight,k,tops,bases)
            %Convert tops and bases to shape matrices.  These are taken as to be the
            %shapes of uniform level/elevation after the rods are pushed.
            
            
            levels = max(tops(:));
            
            for level = 1:levels
                %check if slab is new
                if this.differentSlab(tops,bases,level)
                    k = k+1; %increment number of unique slabs
                    shape_k = (bases <= level).*(level <= tops); %shape of current slab
                    shapes(:,:,k) = shape_k;
                end
                shapesWeight(k) = shapesWeight(k)+1; %if slab is not unique, this increments weight again
            end
        end

        function diffSlab = differentSlab(~,tops,bases,level)

            %Returns 1 if slab level is different than slab level-1 0 otherwise
            
            if level == 1 %first slab is automatically different
                diffSlab = 1;
            else
                shapeLevel = (bases <= level).*(level <= tops); %shape of slab with current level
                shapeLevel_1 = (bases <= level-1).*(level-1 <= tops); %shape of slab with previous level
                diffSlab = ~isequal(shapeLevel,shapeLevel_1); %tests if slabs are equal; isequaln was not giving correct results
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
                checkModality = any(strcmp(matRad_SequencingPhotonsSiochiLeaf.possibleRadiationModes, machine.meta.radiationMode)) && any(strcmp(matRad_SequencingPhotonsSiochiLeaf.possibleRadiationModes, pln.radiationMode));
                
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