function D = matRad_plotDoseCube_4DBio(ct,cst,D, DVH)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots 9 dose cubes based on the dose dsitributions stored in D 
% call  D = matRad_4dRBEcalc(ct, cst, dij, resultGUI, pln, File) to calculate D 
% 
%
% input
%   ct:          matRads ct struct
%   cst:         matRads cst struct
%   D: array with dose cubes
%    e.g.
%   D.name = {'Dopt', 'Drecalc3D', 'Dopt -Drecalc3D', 'Drecalc4Dconst', 'Drecalc4Dvar', 'Drecalc4Dconst-Drecalc4Dvar','Dopt - Drecalc4Dconst', 'Drecalc3D-Drecalc4Dvar', 'Dopt-Drecalc4Dvar'};
%   D.data = {Dopt, Drecalc3D, Dopt-Drecalc3D, Drecalc4Dconst, Drecalc4Dvar, Drecalc4Dconst-Drecalc4Dvar, Dopt - Drecalc4Dconst, Drecalc3D-Drecalc4Dvar, Dopt-Drecalc4Dvar};
%   D.isolines = {1, 1, 0, 1, 1, 0, 0, 0, 0}
%   D.fractions = pln.numOfFractions    
%   DVH: 1 DVH is plotted, 0 no DVH
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
slice = 97; %Liver007
%slice = 110;  %Liver002
vLevels         = [ 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.75 0.8 0.9 0.95 1 1.1 ];


ctSlice      = squeeze(ct.cube{1}(:,:,slice));


%Dosisstatistik
indices     = cst{5,4}{1};  %GTV
indicesPTV = cst{4,4}{1};  %PTV
numOfVoxels = numel(indices);
numOfVoxelsPTV = numel(indicesPTV);

indicesLivermGTV = cst{6,4}{1};  %Liver-GTV

for i = 1:9
    if D.isolines{i} 
    % get Dose, dose is sorted to simplify calculations
    relevantDose = D.data{i};
    doseInVoi    = sort(D.data{i}(indices));
    doseInPTV    = sort(D.data{i}(indicesPTV));
    doseInLiver = sort(D.data{i}(indicesLivermGTV));
    
    refGy = round(linspace(1,round(max(relevantDose(:))),3));
    
    DX = @(x) doseInVoi(ceil((100-x)*0.01*numOfVoxels));
    VX = @(x) numel(doseInVoi(doseInVoi >= x)) / numOfVoxels;
      
    DXPTV = @(x) doseInPTV(ceil((100-x)*0.01*numOfVoxelsPTV));
    VXPTV = @(x) numel(doseInPTV(doseInPTV >= x)) / numOfVoxelsPTV;
    
    D.statD2PTV{i} = DXPTV(2);    
    D.statD95PTV{i} = DXPTV(95);
    D.statD98PTV{i} = DXPTV(98);
    
    D.statD2GTV{i} = DX(2);    
    D.statD95GTV{i} = DX(95);
    D.statD98GTV{i} = DX(98);
        
    D.HIGTV{i} = (DX(2)/DX(98));
    
    D.V95GTV{i} = VX(2*0.95);  %für 2 Gy prescribed
    D.V107GTV{i} = VX(2*1.07);
    
    D.meanLivermGTV{i} = mean(doseInLiver);
    end
end



for i = 1:9
ax{i} = subplot(3,3,i);
ct_rgb = ind2rgb(uint8(63*squeeze(ctSlice)/max(ct.cube{1}(:))),bone);
imagesc(ct_rgb),hold on,
  
doseSlice    = squeeze(D.data{i}(:,:,slice));
minDose         =  min(doseSlice(:));
maxDose         =  max(doseSlice(:));
h1=imagesc(doseSlice);
set(gca,'XTickLabel','');set(gca,'YTickLabel','');
set(h1,'AlphaData', .8*(double(doseSlice)~=0));
%Isolines nicht für Dosisdifferenzen
if D.isolines{i}

vLevelDose      = maxDose  * vLevels;
[~,~] =contour(doseSlice,vLevelDose,'LevelListMode','manual','LineWidth',1.5);  
colormap jet;
else       
Range = linspace(minDose,maxDose,62);
[~,idx]  = min(abs(Range));
idx2 = 62-idx;

a = linspace(0,1,idx);
b = linspace(1,0,idx2);
d1 = ones(1,idx);
d2 = ones(1,idx2);

   blueRow  = [d1 b];
    greenRow = [a b];
    redRow   = [a d2];
    costumMap = [blueRow; greenRow; redRow]'; 
    
colormap(ax{i},costumMap); colorbar;    
end
if D.isolines{i}
%title([D.name{i} 'D95% GTV = '  num2str(D.statD95GTV{i})]) 
title(D.name{i})
else
title(D.name{i})
end


cBarHandel = colorbar(gca);
caxis(gca,[minDose, maxDose]);
%set(get(cBarHandel,'ylabel'),'String', 'dose [Gy]','fontsize','Interpreter','Latex');
set(cBarHandel,'YLim',[minDose maxDose]);

    



colors = colorcube; maskCst = zeros(ct.cubeDim);
    ColorCnt = 1; CntLegend = 1; cHandle = [];
    for i = 1:size(cst,1)
             maskCst(:) = 0;
             maskCst(cst{i,4}{1}) = 1;
            
            if sum(sum(maskCst(:,:,slice))) > 0
                tmpSlice = squeeze(maskCst(:,:,slice));
             else
                tmpSlice = 0;
             end

             if sum(tmpSlice(:)) > 0
                 [c, h] = contour(gca,squeeze(maskCst(:,:,slice)),.5*[1 1],'LineWidth',2,'Color',colors(ColorCnt,:),'DisplayName',cst{i,2});
                 if ~isempty(c)
                     cHandle = [cHandle h]; ColorCnt = ColorCnt + 1; Name = regexprep(cst{i,2},'_','.'); 
                     sLegend{CntLegend} = Name; CntLegend = CntLegend +1;
                 end
             end
    end
    %legend(cHandle,sLegend,'FontSize',defFontSize,'Location','northwestoutside')
end
linkaxes([ax{1}, ax{2}, ax{3}, ax{4}, ax{5}, ax{6}, ax{7}, ax{8}, ax{9}],'xy');



if DVH
figure
n=1000;
c = 1;
legendName{c} = D.name{1};
for i=1:9    
if D.isolines{i} 
dvhPointsOpt = linspace(0,max(D.data{i}(:))*1.05,n);
dvh       = NaN * ones(1,n);
%indices     = cst{5,4}{1};  %GTV
numOfVoxels = numel(indices);
doseInVoi   = D.data{i}(indices); 
for j = 1:n
    dvh(j) = sum(doseInVoi > dvhPointsOpt(j));
end
dvhOpt = dvh ./ numOfVoxels * 100;
  
    plot(dvhPointsOpt, dvhOpt, 'LineWidth', 1.5);hold on
legendName{c} = D.name{i};
    c = c+1;
end
end 
legend(legendName)
title('DVH GTV')
figure
n=1000;
c = 1;
legendName{c} = D.name{1};
for i=1:9    
if D.isolines{i} 
dvhPointsOpt = linspace(0,max(D.data{i}(:))*1.05,n);
dvh       = NaN * ones(1,n);
indicesLiver     = cst{3,4}{1};  %Liver
numOfVoxels = numel(indicesLiver);
doseInVoi   = D.data{i}(indicesLiver); 
for j = 1:n
    dvh(j) = sum(doseInVoi > dvhPointsOpt(j));
end
dvhOpt = dvh ./ numOfVoxels * 100;
  
    plot(dvhPointsOpt, dvhOpt, 'LineWidth', 1.5);hold on
legendName{c} = D.name{i};
    c = c+1;
end
end 
legend(legendName)
title('DVH Liver')
end

figure
n=1000;
c = 1;
legendName{c} = D.name{1};
for i=1:9    
if D.isolines{i} 
dvhPointsOpt = linspace(0,max(D.data{i}(:))*1.05,n);
dvh       = NaN * ones(1,n);
indicesLiver     = cst{6,4}{1};  %Liver -GTV
numOfVoxels = numel(indicesLiver);
doseInVoi   = D.data{i}(indicesLiver); 
for j = 1:n
    dvh(j) = sum(doseInVoi > dvhPointsOpt(j));
end
dvhOpt = dvh ./ numOfVoxels * 100;
  
    plot(dvhPointsOpt, dvhOpt, 'LineWidth', 1.5);hold on
legendName{c} = D.name{i};
    c = c+1;
end
end 
legend(legendName)
title('DVH Liver-GTV')
end

