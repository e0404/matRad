load('baseforweights.mat')

            radius=2;
            ang1=[0,pi/3,pi*2/3,pi,pi*4/3,pi*5/3, ...
                0,pi/6, pi/3,pi/2, pi*2/3,pi*5/6, pi,pi*7/6,...
                pi*4/3,pi*3/2, pi*5/3,pi*11/6];
            nFineSamp = 19;
            for n=1:nFineSamp
                if n==1
                    targetP = stf(i).ray(j).targetPoint_bev;
                elseif n<=7
                    targetP = stf(i).ray(j).targetPoint_bev +...
                        [radius*cos(ang1(n-1)),0,radius*sin(ang1(n-1))];
                else
                    targetP = stf(i).ray(j).targetPoint_bev +...
                        [2*radius*cos(ang1(n-1)),0,2*radius*sin(ang1(n-1))];
                end
                
                [ix,radialDist_sq] = matRad_calcGeoDists(rot_coordsV, ...
                                                     stf(i).sourcePoint_bev, ...
                                                     targetP, ...
                                                     machine.meta.SAD, ...
                                                     radDepthIx, ...
                                                     maxLateralCutoffDoseCalc);
            
            radDepths = radDepthV{1}(ix);   
            
            % just use tissue classes of voxels found by ray tracer
            if (isequal(pln.bioOptimization,'LEMIV_effect') || isequal(pln.bioOptimization,'LEMIV_RBExD')) ... 
                && strcmp(pln.radiationMode,'carbon')
                    vTissueIndex_j = vTissueIndex(ix,:);
            end

                
                % find energy index in base data
                energyIx = find(round2(stf(i).ray(j).energy(k),4) == round2([machine.data.energy],4));
                
                % find depth depended lateral cut off
                if cutOffLevel >= 1
                    currIx = radDepths <= machine.data(energyIx).depths(end) + machine.data(energyIx).offset;
                elseif cutOffLevel < 1 && cutOffLevel > 0
                    % perform rough 2D clipping
                    currIx = radDepths <= machine.data(energyIx).depths(end) + machine.data(energyIx).offset & ...
                         radialDist_sq <= max(machine.data(energyIx).LatCutOff.CutOff.^2);

                    % peform fine 2D clipping  
                    if length(machine.data(energyIx).LatCutOff.CutOff) > 1
                        currIx(currIx) = matRad_interp1((machine.data(energyIx).LatCutOff.depths + machine.data(energyIx).offset)',...
                            (machine.data(energyIx).LatCutOff.CutOff.^2)', radDepths(currIx)) >= radialDist_sq(currIx);
                    end
                else
                    error('cutoff must be a value between 0 and 1')
                end
                
                % empty bixels may happen during recalculation of error
                % scenarios -> skip to next bixel
                if ~any(currIx)
                    continue;
                end
               
                % calculate particle dose for bixel k on ray j of beam i
                %radialDist_sqn = radialDist_sq(:,n);
                if n>1
                tempBixelDose_w{n} = matRad_calcParticleDoseBixel(...
                    radDepths(currIx), ...
                    radialDist_sq(currIx), ...
                    stf(i).ray(j).SSD, ...
                    stf(i).ray(j).focusIx(k), ...
                    machine.data(energyIx));
                [~,idxsIntoTempB] = intersect(superIdx,V(ix(currIx)));
                [~,idxsIntoV] = intersect(V(ix(currIx)),superIdx);
                %disp([size(tempBixelDose,1) max(idxsIntoV) size(bixelDose,1) max(idxsIntoTempB)]);
                
                bixelDose_w(idxsIntoTempB) = bixelDose_w(idxsIntoTempB) + tempBixelDose_w{n}(idxsIntoV);
%                 idc = V(ix(currIx)); idc(idxsIntoV)=[];
%                 totidx = cat(1,tempBixelIdx,idc);
%                 [~,idxsIntoC] = intersect(V(ix(currIx)),idc);
%                 bixelDose = cat(1,bixelDose,tempBixelDose(idxsIntoC));
%                 [superIdx,sortidx]=sort(totidx);
%                 bixelDose = bixelDose(sortidx);
                else
                    bixelDose_w = matRad_calcParticleDoseBixel(...
                    radDepths(currIx), ...
                    radialDist_sq(currIx), ...
                    stf(i).ray(j).SSD, ...
                    stf(i).ray(j).focusIx(k), ...
                    machine.data(energyIx));
                    superIdx = V(ix(currIx));
                    tempBixelDose_w{1} = bixelDose_w;
                end
                tempBixelIdx = V(ix(currIx));
                doseTmpContainer_w{n}= sparse(V(ix(currIx)),1,tempBixelDose_w{n},dij.numOfVoxels,1);
                doseTmpContainer_w{n}= full(doseTmpContainer_w{n});
                %doseTmpContainer_w{n}=reshape(doseTmpContainer_w{n},[160 160 160]);
            end
            
           
%%
load('baseforweights3.mat');
                doseTmpContainer{1}= full(doseTmpContainer{1});
                %doseTmpContainer{1}=reshape(doseTmpContainer{1},[160 160 160]);


for i=1:size(superIdx,1);
    totBixelDose = @(x) doseTmpContainer{1}(superIdx(i)) - ( x(1).*( doseTmpContainer_w{1}(superIdx(i)) + ...
        exp(-radius/(2*x(2))).*( doseTmpContainer_w{2}(superIdx(i)) + doseTmpContainer_w{3}(superIdx(i)) +...
        doseTmpContainer_w{4}(superIdx(i)) + doseTmpContainer_w{5}(superIdx(i)) + doseTmpContainer_w{6}(superIdx(i)) +...
        doseTmpContainer_w{7}(superIdx(i)) ) + exp(-(2*radius)/(2*x(2))).*...
        ( doseTmpContainer_w{8}(superIdx(i)) + doseTmpContainer_w{9}(superIdx(i)) + doseTmpContainer_w{10}(superIdx(i)) +...
        doseTmpContainer_w{11}(superIdx(i)) + doseTmpContainer_w{12}(superIdx(i)) + doseTmpContainer_w{13}(superIdx(i)) +...
        doseTmpContainer_w{14}(superIdx(i)) + doseTmpContainer_w{15}(superIdx(i)) + doseTmpContainer_w{16}(superIdx(i)) +...
        doseTmpContainer_w{17}(superIdx(i)) + doseTmpContainer_w{18}(superIdx(i)) + doseTmpContainer_w{19}(superIdx(i))) ) );
    
    [xfin(:,i),fval] = fminsearch(totBixelDose,[.1, .1]);
    if mod(i,1000)==0
        disp(i)
    end
end

    
    
%     
%     
%     