function h = plot_DCH(voiName,voiVolume,cst,doseVec,dij,color)

LineWidth = 1.5;
LineStyle = {'-','--',':','-.'};

% plot DCH
for i = 1:length(voiName)
    % plot DCH
    [dchPoints,Q] = matRad_calcDCH(voiVolume(i),doseVec,dij,cst(strcmp([cst(:,2)],voiName{i}),:));
    h(i,1)        = plot(dchPoints,Q,LineStyle{i},'Color',color,'LineWidth',LineWidth);
    hold on
    
    % plot references if available
    cstidxVOIScenUnion = find(strcmp([cst(:,2)],[voiName{i},' ScenUnion']));
   
    if ~isempty(cstidxVOIScenUnion)
        if ~isempty(cst{cstidxVOIScenUnion,6})
            logidxDCHobj = ~cellfun('isempty',strfind({cst{cstidxVOIScenUnion,6}(:).type},'DCH'));
            volume       = [cst{cstidxVOIScenUnion,6}(logidxDCHobj).volume]./100;
            coverage     = [cst{cstidxVOIScenUnion,6}(logidxDCHobj).coverage]./100;
            dose         = [cst{cstidxVOIScenUnion,6}(logidxDCHobj).dose];
            for j = 1:length(coverage)
                if isequal(voiVolume(i),volume(j))
                    h(i,2)       = scatter(dose(j),coverage(j)*100,60,'o','MarkerEdgeColor','k','LineWidth',LineWidth);
                end
            end
        end
    end
    
    cstidxVOI = find(strcmp([cst(:,2)],voiName{i}));
    
    if ~isempty(cstidxVOI)
        if ~isempty(cst{cstidxVOI,6})
            logidxDCHobj = ~cellfun('isempty',strfind({cst{cstidxVOI,6}(:).type},'DCH'));
            volume       = [cst{cstidxVOI,6}(logidxDCHobj).volume]./100;
            coverage     = [cst{cstidxVOI,6}(logidxDCHobj).coverage]./100;
            dose         = [cst{cstidxVOI,6}(logidxDCHobj).dose];
            for j = 1:length(coverage)
                if isequal(voiVolume(i),volume(j))
                    h(i,2)       = scatter(dose(j),coverage(j)*100,60,'o','MarkerEdgeColor','k','LineWidth',LineWidth);
                end
            end
        end
    end
    
end

% set plot options
xlabel('dose [Gy]')
ylabel('coverage probability [%]')
set(gca,'ygrid','on')
ylim([0 110])
axis square

end