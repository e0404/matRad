function [dchPoints,Q] = matRad_calcDCH(volume,doseVec,cst,numOfScenarios)

n         = 10000;
numOfVois = size(cst,1);

% set DCH points
dchPoints = linspace(0,max(vertcat(doseVec{:}))*1.05,n);

%figure
for i = 1:numOfVois
    indices   = cst{i,4}{1};

    % calculate dose that corresponds to volume and deviation from D
    for Scen = 1:numOfScenarios
        dose = matRad_calcInversDVH(volume(i)/100,doseVec{Scen}(indices));
        dev(Scen,:) = dose - dchPoints;
    end

    % calculate logical mask
    devlog = dev >= 0;

    % calculate coverage Q
    Q(i,:) = (1/numOfScenarios)*sum(devlog)*100;

    % plot coverage over dose
    %plot(dchPoints,Q(i,:))
    %hold on

    % store legend info
    %legendinfo{i} = [cst{i,2},' D_{',num2str(volume(i)),'} DCH'];
end

% set legend
%legend(legendinfo)
%grid on
%xlabel('Dose [Gy]')
%ylabel('Coverage Probability [%]')

end