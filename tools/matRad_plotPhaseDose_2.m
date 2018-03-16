function matRad_plotPhaseDose_2(resultGUI, slice)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Silke
% plots dose cubes for each phase
% 
% call
%   matRad_plotPhaseDose(slice)
%
% input
%   
% output
%
%
% References
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% get data from workspace
cst     = evalin('base','cst');
ct     = evalin('base','ct');
pln    = evalin('base','pln');

if(nargin < 2)
% get isocenter slices
vIsoCenter = round(pln.isoCenter./[ct.resolution.x ct.resolution.y ct.resolution.z]);
slice = vIsoCenter(3);
end


%%
% number of plots
NumberPlots = length(resultGUI.phaseDose);
nr = ceil(sqrt(NumberPlots));
nc = ceil(NumberPlots / nr);

figure
set(gcf,'Color',[1 1 1]);
n = 1;
    
%Dose in each phase
while( n < length(resultGUI.phaseDose)+1)
subplot(nr,nc,n)
imagesc(resultGUI.phaseDose{n}(:,:,slice))
colorbar
title(['Dose in CT phase ', num2str(n)])
n= n+1;
end

%transformed dose
% figure
% set(gcf,'Color',[1 1 1]);
% n = 1;
% while( n < length(resultGUI.phaseDose)+1)
% % transversal slices
% subplot(nr,nc,n)
% imagesc(resultGUI.phaseDoseDDM{n}(:,:,slice))
% colorbar
% title(['transf dose in CT phase ', num2str(n)])
% n= n+1;
% end




figure
subplot(3,2,1)

 if  strcmp(pln.bioOptimization,'none')
imagesc(resultGUI.physicalDose(:,:,slice));
Int = sum(resultGUI.physicalDose(:));
 else
imagesc(resultGUI.RBExD(:,:,slice));
Int = sum(resultGUI.RBExD(:));
 end

colorbar
title(['optimized dose in slice' num2str(slice) ', Int Dose = ' num2str(Int)])
    
%accumulated dose
subplot(3,2,3)
 if  strcmp(pln.bioOptimization,'none')
imagesc(resultGUI.accDose(:,:,slice));
Int = sum(resultGUIaccDose(:));
 else
imagesc(resultGUI.accRBExD(:,:,slice));
Int = sum(resultGUI.accRBExD(:));
 end
 
colorbar
title(['accumulated dose Intdose = ' num2str(Int)])

% %diff dose
subplot(3,2,5)
 if  strcmp(pln.bioOptimization,'none')
imagesc(resultGUI.physicalDose(:,:,slice) - resultGUI.accDose(:,:,slice));
 else
imagesc(resultGUI.RBExD(:,:,slice) - resultGUI.accRBExD(:,:,slice));
 end 
colorbar
title('opt - acc ')


%calc DVH final dose
n=1000;
 if  strcmp(pln.bioOptimization,'none')
dvhPoints = linspace(0,max(resultGUI.physicalDose(:))*1.05,n);
 else
dvhPoints = linspace(0,max(resultGUI.RBExD(:))*1.05,n);
 end 
 %dvhPoints = linspace(0,max(resultGUI.finalDose(:))*1.05,n);
dvh       = NaN * ones(1,n);
subplot(3,2,2)
numOfVois = size(cst,1);
for i = 1:numOfVois
    if cst{i,5}.Visible
        indices     = cst{i,4}{1};
        numOfVoxels = numel(indices);
         if  strcmp(pln.bioOptimization,'none')
 doseInVoi   = resultGUI.physicalDose(indices);
 else
doseInVoi   = resultGUI.RBExD(indices);
 end 
       

        for j = 1:n
            dvh(j) = sum(doseInVoi > dvhPoints(j));
        end

        dvh = dvh ./ numOfVoxels * 100;
        plot(dvhPoints, dvh);hold on
    end
end

%calc DVH acc dose
 if  strcmp(pln.bioOptimization,'none')
dvhPoints = linspace(0,max(resultGUI.accDose(:))*1.05,n);
 else
dvhPoints = linspace(0,max(resultGUI.accRBExD(:))*1.05,n);
 end 
dvh       = NaN * ones(1,n);

subplot(3,2,4)
numOfVois = size(cst,1);
for i = 1:numOfVois
    if cst{i,5}.Visible
        indices     = cst{i,4}{1};
        numOfVoxels = numel(indices);
                if  strcmp(pln.bioOptimization,'none')
 doseInVoi   = resultGUI.accDose(indices);
 else
doseInVoi   = resultGUI.accRBExD(indices);
 end 

        for j = 1:n
            dvh(j) = sum(doseInVoi > dvhPoints(j));
        end

        dvh = dvh ./ numOfVoxels * 100;
        plot(dvhPoints, dvh);hold on
    end
end


 


