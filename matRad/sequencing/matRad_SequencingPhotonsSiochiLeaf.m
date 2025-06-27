classdef  matRad_SequencingPhotonsSiochiLeaf < matRad_SequencingPhotonsAbstract

    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties (Constant)
        name = 'Photons Sochi Leaf Sequenceer';
        shortName = 'Sochi Leaf';
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
                sequence.beam(i).shapesWeight = shapesWeight(1:k)/this.numOfLevels*calFac;
                sequence.beam(i).bixelIx      = 1+offset:stf(i).numOfRays+offset;
                sequence.beam(i).fluence      = D_0;
                sequence.beam(i).sum          = zeros(size(D_0));
                
                for j = 1:k
                    sequence.beam(i).sum = sequence.beam(i).sum+sequence.beam(i).shapes(:,:,j)*sequence.beam(i).shapesWeight(j);
                end
                sequence.w(1+offset:stf(i).numOfRays+offset,1) = sequence.beam(i).sum(indInMx);
                
                offset = offset + stf(i).numOfRays;
                
            end

            if this.visBool
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

        function plotSegments(this,sequencing)
                % create the sequencing figure
                sz = [800 1000]; % figure size
                screensize = get(0,'ScreenSize');
                xpos = ceil((screensize(3)-sz(2))/2); % center the figure on the screen horizontally
                ypos = ceil((screensize(4)-sz(1))/2); % center the figure on the screen vertically
                seqFig = figure('position',[xpos,ypos,sz(2),sz(1)]);  

                for i = 1:numel(sequencing)

                    D_0 = sequencing.beam(i).fluence;

                    clf(seqFig);
                    colormap(seqFig,'jet');
                        
                    seqSubPlots(1) = subplot(2,2,1,'parent',seqFig);
                    imagesc(sequencing.beam(i).fluence,'parent',seqSubPlots(1));
                    set(seqSubPlots(1),'CLim',[0 this.numOfLevels],'YDir','normal');
                    title(seqSubPlots(1),['Beam # ' num2str(i) ': max(D_0) = ' num2str(max(D_0(:))) ' - ' num2str(numel(unique(D_0))) ' intensity levels']);
                    xlabel(seqSubPlots(1),'x - direction parallel to leaf motion ')
                    ylabel(seqSubPlots(1),'z - direction perpendicular to leaf motion ')
                    colorbar;
                    drawnow

                    %show the leaf positions
                    D_k =  sequencing.beam(i).fluence;
                    for  k = 1:sequencing.beam(i).numOfShapes
                        shape_k = sequencing.beam(i).shapes(:,:,k);
                        [dimZ,dimX] = size(sequencing.beam(i).fluence);
                        seqSubPlots(4) = subplot(2,2,3.5,'parent',seqFig);
                        imagesc(shape_k,'parent',seqSubPlots(4));
                        hold(seqSubPlots(4),'on');
                        set(seqSubPlots(4),'YDir','normal')
                        xlabel(seqSubPlots(4),'x - direction parallel to leaf motion ')
                        ylabel(seqSubPlots(4),'z - direction perpendicular to leaf motion ')
                        title(seqSubPlots(4),['beam # ' num2str(i) ' shape # ' num2str(k) ' d_k = ' num2str(sequencing.beam(i).shapesWeight(k))]);
                        for j = 1:dimZ
                            leftLeafIx = find(shape_k(j,:)>0,1,'first');
                            rightLeafIx = find(shape_k(j,:)>0,1,'last');
                            if leftLeafIx > 1
                                plot(seqSubPlots(4),[.5 leftLeafIx-.5],j-[.5 .5] ,'w','LineWidth',2)
                                plot(seqSubPlots(4),[.5 leftLeafIx-.5],j+[.5 .5] ,'w','LineWidth',2)
                                plot(seqSubPlots(4),[ leftLeafIx-.5 leftLeafIx-.5],j+[.5 -.5] ,'w','LineWidth',2)
                            end
                            if rightLeafIx<dimX
                                plot(seqSubPlots(4),[dimX+.5 rightLeafIx+.5],j-[.5 .5] ,'w','LineWidth',2)
                                plot(seqSubPlots(4),[dimX+.5 rightLeafIx+.5],j+[.5 .5] ,'w','LineWidth',2)
                                plot(seqSubPlots(4),[ rightLeafIx+.5 rightLeafIx+.5],j+[.5 -.5] ,'w','LineWidth',2)
                            end
                            if isempty(rightLeafIx) && isempty (leftLeafIx)
                                plot(seqSubPlots(4),[dimX+.5 .5],j-[.5 .5] ,'w','LineWidth',2)
                                plot(seqSubPlots(4),[dimX+.5 .5],j+[.5 .5] ,'w','LineWidth',2)
                                plot(seqSubPlots(4),.5*dimX*[1 1]+[0.5],j+[.5 -.5] ,'w','LineWidth',2)
                            end
                        end
                        pause(1);
                        
                        %Plot residual intensity matrix.
                        D_k = D_k-shape_k; %residual intensity matrix for visualization
                        seqSubPlots(2) = subplot(2,2,2,'parent',seqFig);
                        imagesc(D_k,'parent',seqSubPlots(2));
                        set(seqSubPlots(2),'CLim',[0 this.numOfLevels],'YDir','normal');
                        title(seqSubPlots(2),['k = ' num2str(k)]);
                        colorbar
                        drawnow
                        
                        axis tight
                        drawnow
                    end


                end
        end

    end
end