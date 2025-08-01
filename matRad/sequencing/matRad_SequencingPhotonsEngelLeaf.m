classdef  matRad_SequencingPhotonsEngelLeaf < matRad_SequencingPhotonsAbstract

% multileaf collimator leaf sequencing algorithm 
% for intensity modulated beams with multiple static segments accroding 
% to Engel et al. 2005 Discrete Applied Mathematics
%
% References
%   [1] http://www.sciencedirect.com/science/article/pii/S0166218X05001411
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
        name = 'Photons Engel Leaf Sequenceer';
        shortName = 'engel';
        possibleRadiationModes = {'photons'};
    end 

    methods

        function sequence = sequence(this,w,stf)
            matRad_cfg = MatRad_Config.instance();
            numOfBeams = numel(stf);
            offset = 0;

            for i = 1:numOfBeams

                [D_0,D_k, shapes,calFac,indInMx] = this.initBeam(stf(i),w(1+offset:stf(i).numOfRays+offset));

                 k = 0;
                
                % start sequencer
                while max(D_k(:) > 0)
                    
                    %calculate the difference matrix diffMat
                    diffMat = diff([zeros(size(D_k,1),1) D_k zeros(size(D_k,1),1)],[],2);
                    
                    %calculate complexities
                    c = sum(max(0,diffMat),2); %TNMU-row-complexity
                    com = max(c); %TNMU complexity
                    g = com - c; %row complexity gap
                    
                    %initialize segment
                    segment = zeros(size(D_k));
                    
                    k = k + 1;
                                       
                    
                    %loop over all rows 
                    for j=1:size(D_0,1)
                        
                        %determine essential intervals
                        data(j).left(1) = 0; %left interval limit, actual for an empty interval
                        data(j).right(1) = 0; %right interal limit, actual for an empty interval
                        data(j).v(1) = g(j);  %greatest number such that the inequalities (6) resp. (7) is satisfied with u=v
                        data(j).w(1) = inf; %smallest number in the interval
                        data(j).u(1) = data(j).v(1); %min(v,w)
                        
                        [~, pos, ~] = find(diffMat(j,:) > 0); % indices of all positive elements in the j. row of diffmat
                        [~, neg, ~] = find(diffMat(j,:) < 0); % indices of all negative elements in the j. row of diffMat
                        
                        n=2;
                        
                        %loop over the positive elements in the j. row of diffmat ->
                        %possible left interval limits
                        for m=1:size(pos,2)
                            
                            %loop over the negative elements in the j. row of diffMat ->
                            %possible right interval limit
                            for l=1:size(neg,2)
                                
                                %take only intervals I=[l,r] with l<=r
                                if pos(m) <= neg(l)-1
                                    
                                    %set interval limits
                                    data(j).left(n) = pos(m);
                                    data(j).right(n) = neg(l)-1;
                                    
                                    %calculate v according to Lemma 8
                                    if g(j) <= abs( diffMat(j,pos(m)) + diffMat(j,neg(l)) )
                                        data(j).v(n) = min( diffMat(j,pos(m)), -diffMat(j,neg(l)) ) + g(j);
                                    else
                                        data(j).v(n) = ( diffMat(j, pos(m)) - diffMat(j, neg(l)) + g(j)) / 2;
                                    end
                                    
                                    %calculate w and u according to equality (11) and
                                    %(12) 
                                    data(j).w(n) = min(D_k(j,pos(m):(neg(l)-1)));
                                    data(j).u(n) = min(data(j).v(n), data(j).w(n));
                                    
                                    n = n+1;
                                end
                            end
                        end
                        
                        u(j) = max(data(j).u);
                        
                    end
                    
                    %calculate u_max from theorem 9
                    d_k = min(u);
                    
                    %loop over all rows
                    for j=1:size(D_0,1)
             
                        %find all possible (and essential) intervals
                        candidate = find(data(j).u >= d_k); 
                        
                        %calculate the potential of the possible intervals
                        
                        %initialize p as -Inf
                        data(j).p(1:length(data(j).left)) = -Inf;
                        
                        %loop over all possible intervals
                        for s=1:size(candidate,2)
                            
                            if (s==1 && data(j).left(candidate(s)) == 0)
                                data(j).p(candidate(1)) = 0;
                                
                                
                            else
                                %calculate p1 according to equality (17)
                                if (d_k == diffMat(j, data(j).left(candidate(s))) && d_k ~= D_k(j, data(j).left(candidate(s))))
                                    p1 = 1;
                           
                                else
                                    p1 = 0;
                                
                                end
                                
                                %calculate p2 according to equalitiy (18)
                               % if data(j).right(candidate(s)) < size(D_0, 2)
                                    
                                    if (d_k == -diffMat(j, data(j).right(candidate(s))+1) && d_k ~= D_k(j, data(j).right(candidate(s))))
                                        p2 = 1;
                                    else
                                        p2 = 0;
                                    end
                                    
            %                     else
            %                         
            %                         if d_k == -diffMat(j, data(j).right(candidate(s))+1)
            %                             p2 = 1;
            %                         else
            %                             p2 = 0;
            %                         end
            %                         
            %                     end                                     
                                
                                %calculate p3 according to equality (19)
                                p3 = size(find(D_k(j, data(j).left(candidate(s)):data(j).right(candidate(s))) == d_k),2);
                                
                                data(j).p(candidate(s)) = p1 + p2+ p3;
                                
                            end
                            
                        end
                            
                            %determinate intervals with maximum potential
                            maxPot = find(data(j).p == max(data(j).p));
                            
                            %if several intervals have maximum potential, select
                            %the interval which has maximum length
                            if size(maxPot,2) > 1
            
                                for t=1:size(maxPot,2)
                                    if t==1 && data(j).left(maxPot(t)) == 0
                                        data(j).l(1) = 0;
                                    else
                                        data(j).l(maxPot(t)) = data(j).right(maxPot(t)) - data(j).left(maxPot(t)) + 1;
                                    end
                                end
                                
                                %data(j).l(maxPot) = data(j).right(maxPot) - data(j).left(maxPot) + 1;
                                
                                maxLength = find(data(j).l == max(data(j).l));
                                 
                                %left and right interval limits of the selected
                                %interval
                                leftIntLimit(j) = data(j).left(maxLength(1));
                                rightIntLimit(j) = data(j).right(maxLength(1));
                                
             
                            else
                                
                                %left and right interval limits of the selected
                                %interval
                                leftIntLimit(j) = data(j).left(maxPot);
                                rightIntLimit(j) = data(j).right(maxPot);
                                
                                
                            end
                            
                            %create segment associated by the selected interval
                            if leftIntLimit(j) ~= 0 
                                
                                segment(j,leftIntLimit(j):rightIntLimit(j)) = 1;
             
                            end
            
                    end
                    
                    %write the segment in shape_k
                    shape_k = segment;                                    
            
                    %save shape_k in container
                    shapes(:,:,k) =  shape_k;
                    
                    %save the calculated MU
                    shapesWeight(k) = d_k; 
                    
                    %calculate  new matrix, the  diference matrix and complexities
                    D_k = D_k - d_k*shape_k; 
                    
                    %delete variables
                    clear data;
                    clear segment;
                    clear u;
                    clear leftIntLimit;
                    clear rightIntLimit;
                   
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
                checkModality = any(strcmp(matRad_SequencingPhotonsEngelLeaf.possibleRadiationModes, machine.meta.radiationMode)) && any(strcmp(matRad_SequencingPhotonsEngelLeaf.possibleRadiationModes, pln.radiationMode));
                
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