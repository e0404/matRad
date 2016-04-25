function matRad_calcDCH(volume,doseVec,cst,numOfScenarios)
%[Q,D] = matRad_calcDCH(volume,doseVec,cst,multScen)

numOfVois = size(cst,1);

% set dose points D
D = linspace(0,max(vertcat(doseVec{:}))*1.05,10000);

figure
for i = 1:numOfVois
    if ~isequal(cst{i,2},'TargetRing')
        indices   = cst{i,4}{1};

        % calculate dose that corresponds to volume and deviation from D
        for Scen = 1:numOfScenarios
            dose = matRad_calcInversDVH(volume(i),doseVec{Scen}(indices));
            dev(Scen,:) = dose - D;
        end

        % calculate logical mask
        devlog = dev >= 0;

        % calculate coverage Q
        Q(i,:) = (1/numOfScenarios)*sum(devlog);

        % plot coverage over dose
        plot(D,Q(i,:))
        hold on

        % store legend info
        legendinfo{i} = [cst{i,2},' D_{',num2str(volume(i)*100),'} DCH'];
    end
end

% set legend
legend(legendinfo)
grid on

end