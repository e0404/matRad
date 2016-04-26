function matRad_calcPDVH(Coverage,doseVec,cst,numOfScenarios)

n         = 10000;
numOfVois = size(cst,1);

% set DVH points
dvhPoints = linspace(0,max(vertcat(doseVec{:}))*1.05,n);

figure
for i = 1:numOfVois
    dvh = repmat(NaN * ones(1,n),numOfScenarios,1);
    for Scen = 1:numOfScenarios
        indices     = cst{i,4}{1};
        numOfVoxels = numel(indices);
        doseInVoi   = doseVec{Scen}(indices);

        % calculate DVH
        for j = 1:n
            dvh(Scen,j) = sum(doseInVoi > dvhPoints(j));
        end

        dvh(Scen,:) = dvh(Scen,:) ./ numOfVoxels * 100;
    end

    % calculate PDVH
    for j = 1:n
        VolumePointsSorted = sort(dvh(:,j),'descend');
        ix = max([1 ceil(Coverage/100*numel(VolumePointsSorted))]);
        V(j) = VolumePointsSorted(ix);
    end
    
    plot(dvhPoints,V)
    
    % store legend info
    legendinfo{i} = [cst{i,2},' PDVH_{',num2str(Coverage(i)),'}'];

end

% set legend
legend(legendinfo)
grid on
xlabel('Dose [Gy]')
ylabel('Volume [%]')

end