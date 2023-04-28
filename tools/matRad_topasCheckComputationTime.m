function computationTimeInHours = matRad_topasCheckComputationTime(folder)

d = dir(folder);
d = d(~[d.isdir]);

startedProcessing = contains({d.name},'out','IgnoreCase',true) & ~contains({d.name},'_std','IgnoreCase',true);
startDate = [d(startedProcessing).datenum];

finishedProcessing = contains({d.name},'physicalDose.binheader','IgnoreCase',true) & ~contains({d.name},'_std','IgnoreCase',true);
finishDate = [d(finishedProcessing).datenum];

% Calculate difference and convert to hours
computationTimeInHours = 24 * sum(finishDate - startDate);

hours = floor(computationTimeInHours);
minutes = round((computationTimeInHours - hours)*60);

disp(['Simulation time was ' num2str(hours), 'h:' num2str(minutes) 'm with ' num2str(sum(finishedProcessing)) ' batches.'])

end

