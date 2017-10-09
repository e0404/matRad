function matRad_latexReport(ct, cst, pln, nominalScenario, structureStat, doseStat, resultGUI, param)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad uncertainty analysis report generaator function
% 
% call
%   latexReport(ct, cst, pln, nominalScenario, structureStat, param)
%
% input
%   ct:                 ct cube
%   cst:                matRad cst struct
%   pln:                matRad plan meta information struct
%   nominalScenario:    struct containing dose, qi and dvh of the nominal scenario
%   structureStat:      structures which were examined (can be empty, 
%                       when all structures were examined)
%   param               (optional) struct set of parameters, such as output path

% output
%   (binary)            a pdf report will be generated and saved
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y, argmin] = cutAtArgmin(x)
  [~,argmin] = min(x);
  y = x(1:argmin);
end

if exist('param','var') && ~isempty(param)
    outputPath = param.reportPath;
end


%% create latex Report
if ~exist('outputPath','var')
    outputPath = uigetdir;
end
%outputPath = fullfile('tools','reportGeneration','latex','output');


%% correct cst for unwanted characters and disable commonly not wanted structures
notVisibleStructs = {'Beekleys', 'Beekley', 'CT-Referenzpunkt'};
for i = 1:size(cst,1)
    % use only alphanumerical characters
    cst{i,2} = regexprep(cst{i,2},'[^a-zA-Z0-9]','-');
    if isempty(cst{i,4}{1}) || (sum(strcmp(cst{i,2}, notVisibleStructs)) >= 1)
        cst{i,5}.Visible = false;
    end
end


%% insert standard patient information
try
    % import patient information
    patientInformation.firstName = ct.dicomInfo.PatientName.GivenName;
    patientInformation.lastName = ct.dicomInfo.PatientName.FamilyName;
    patientInformation.sex = ct.dicomMeta.PatientSex;
    patientInformation.patientID = ct.dicomMeta.PatientID;
catch
    patientInformation.firstName = 'N.A.';
    patientInformation.lastName = 'N.A.';
    patientInformation.sex = 'N.A.';
    patientInformation.patientID = 'N.A.';
end

% import plan information
planInformation.gantryAngles = num2str(pln.gantryAngles);
planInformation.couchAngles = num2str(pln.couchAngles);
planInformation.modality = pln.radiationMode;

line = cell(0);
line =  [line; '\newcommand{\patientFirstName}{',patientInformation.firstName,'}'];
line =  [line; '\newcommand{\patientLastName}{',patientInformation.lastName,'}'];
line =  [line; '\newcommand{\patientSex}{',patientInformation.sex,'}'];
line =  [line; '\newcommand{\patientID}{',patientInformation.patientID,'}'];
line =  [line; '\newcommand{\operator}{',param.operator,'}'];

line =  [line; '\newcommand{\reportGenerationDate}{\today}'];

line =  [line; '\newcommand{\planGantryAngles}{',planInformation.gantryAngles,'}'];
line =  [line; '\newcommand{\planCouchAngles}{',planInformation.couchAngles,'}'];
line =  [line; '\newcommand{\planRadiationModality}{',planInformation.modality,'}'];

fid = fopen(fullfile(outputPath,'patientInformation.tex'),'w');
for i = 1:numel(line)
    text = regexprep(line{i},'\','\\\');
    fprintf(fid,text);
    fprintf(fid,'\n');
end
fclose(fid);


%% plot isocentre slices (nominal, mean, std)
for plane=1:3
    switch plane 
        case 1
            slice = round(pln.isoCenter(1,plane) / ct.resolution.x,0);
        case 2
            slice = round(pln.isoCenter(1,plane) / ct.resolution.y,0);
        case 3
            slice = round(pln.isoCenter(1,plane) / ct.resolution.z,0);
    end
    colors = colorcube(size(cst,1));
    for cubesToPlot = 1:3
        if cubesToPlot == 1
            if isfield(resultGUI,'RBExDose')
                doseCube = resultGUI.RBExDose;
                colorMapLabel = 'physical Dose [Gy]';
            else
                doseCube = resultGUI.physicalDose;
                colorMapLabel = 'RBExDose [Gy(RBE)]';
            end
            fileSuffix = 'nominal';
        elseif cubesToPlot == 2
            doseCube = doseStat.meanCubeW;
            fileSuffix = 'meanW';
        elseif cubesToPlot == 3
            doseCube = doseStat.stdCubeW;
            fileSuffix = 'stdW';
        end
        figure; ax = gca;
        matRad_plotSliceWrapper(ax,ct,cst,1,doseCube,plane,slice,[],[],colors,[],colorMapLabel);
        drawnow();
        cleanfigure();          
        matlab2tikz(fullfile(outputPath,['isoSlicePlane', num2str(plane), '_', fileSuffix, '.tex']), 'relativeDataPath', 'data', 'showInfo', false, 'width', '\figW');
        close
    end
end


%% plot nominal dvh and write qi table
% plot dvh
colors = jet(size(cst,1));
hold off;
for i = 1:size(cst,1)
    if cst{i,5}.Visible == true
        [y, argmin] = cutAtArgmin(nominalScenario.dvh{i}(2,:));
        x = nominalScenario.dvh{i}(1,1:argmin);        
        h(1) = plot(x,y,'LineWidth',2, 'Color', colors(i,:), 'DisplayName', cst{i,2});      
        ylim([0 100]);
        if strcmp(pln.bioParam, 'RBExDose')
            xlabel('Dose RBE x [Gy]');
        else
            xlabel('Dose [Gy]');
        end
        ylabel('Volume [%]');
        lh = legend('show','Location','northeastoutside');
        hold on;
    end
end
drawnow;
matlab2tikz(fullfile(outputPath,'nominalDVH.tex'),'showInfo', false, 'width', '0.7\textwidth');
hold off
close

% write qi table

% get list of all fields
fieldsInQi = fieldnames(nominalScenario.qi{1});
for i= 1:numel(fieldsInQi)
    nanRow.(fieldsInQi{i}) = NaN;
end

nomQi = nominalScenario.qi{1};
structName = {};
c = 1;
for i = 1:size(cst,1)
    if cst{i,5}.Visible == true && ~isempty(cst{i,4}{1})
        nomQi(c) = nominalScenario.qi{i};
        fullQi(i) = nominalScenario.qi{i};
        structName = [structName, cst{i,2}];
        c = c + 1;
    else
        fullQi(i) = nanRow;
    end
end

nomQiTable = struct2table(nomQi);
fullNomQiTable = struct2table(fullQi);
fullNomQiTable = fullNomQiTable(:,1:8);
nomQiTable.Properties.RowNames = structName;
input.data = nomQiTable(:,1:8);
input.dataFormat = {'%.2f'};
input.booktabs = 1;
input.tableBorders = 0;
input.tableEnvironment = 0;

latex = latexTable(input);

% save LaTex code as file
filename{i}.QI = regexprep([cst{i,2},'_QI.tex'], '\s+', '');
fid=fopen(fullfile(outputPath,'nominalQI.tex'),'w');
[nrows,ncols] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fclose(fid);
clear latex
clear input


%% analysis parameters
line = cell(0);
if pln.multScen.numOfCtScen <= 1
    line =  [line; '\newcommand{\ctScen}{false}'];
else
    line =  [line; '\newcommand{\ctScen}{true}'];
end

if pln.multScen.numOfRangeShiftScen <= 1
    line =  [line; '\newcommand{\rangeScen}{false}'];
    line =  [line; '\newcommand{\rangeRelSD}{', num2str(0), '}'];
    line =  [line; '\newcommand{\rangeAbsSD}{', num2str(0), '}'];
else
    line =  [line; '\newcommand{\rangeScen}{true}'];
    if isempty(pln.multScen.relRangeShift) || ((max(pln.multScen.relRangeShift) == 0) && (min(pln.multScen.relRangeShift) == 0))
        line =  [line; '\newcommand{\rangeRelSD}{', num2str(0), '}'];
    else
        line =  [line; '\newcommand{\rangeRelSD}{', num2str(pln.multScen.rangeRelSD), '}'];
    end
    if isempty(pln.multScen.relRangeShift) || ((max(pln.multScen.absRangeShift) == 0 && min(pln.multScen.absRangeShift) == 0))
        line =  [line; '\newcommand{\rangeAbsSD}{', num2str(0), '}'];
    else
        line =  [line; '\newcommand{\rangeAbsSD}{', num2str(pln.multScen.rangeAbsSD), '}'];
    end
end

if pln.multScen.numOfShiftScen <= 1
    line =  [line; '\newcommand{\shiftScen}{false}'];
    line =  [line; '\newcommand{\shiftSD}{', num2str([0 0 0]), '}'];
else
    line =  [line; '\newcommand{\shiftScen}{true}'];
    line =  [line; '\newcommand{\shiftSD}{', num2str(pln.multScen.shiftSD), '}'];
end

fid = fopen(fullfile(outputPath,'uncertaintyParameters.tex'),'w');
for i = 1:numel(line)
    text = regexprep(line{i},'\','\\\');
    fprintf(fid,text);
    fprintf(fid,'\n');
end
fclose(fid);


%% add DVH and QI
clear h
clear filename
% relative file path (relative to main.tex)
relativePath = fullfile('data','structures');

for i = 1:size(cst,1)
    if cst{i,5}.Visible == true
        numOfConf = floor(size(structureStat(i).dvhStat.percDVH,1) / 2);
        % DVH
        doseGrid = structureStat(i).dvhStat.mean(1,:);

        % plot nominal plan
        [y, argmin] = cutAtArgmin(nominalScenario.dvh{i}(2,:));
        x = nominalScenario.dvh{i}(1,1:argmin); 
        h(1) = plot(x,y,'LineWidth',2, 'Color', 'k', 'DisplayName', 'nominal');
        hold on;
        % plot mean
        [y, argmin] = cutAtArgmin(structureStat(i).dvhStat.mean(2,:));
        x = structureStat(i).dvhStat.mean(1,1:argmin); 
        h(2) = plot(x,y,'--','LineWidth',2, 'Color', 'k', 'DisplayName', '\mu');
        % plot dvh confidence bands
        % colors
        colors = jet(numOfConf);
        alphaTrans = 1;

        hIx = numel(h);
        for j = 1:numOfConf
            hIx = hIx + 1;
            lIx = j;
            hIx = size(structureStat(i).dvhStat.percDVH,1) - (j-1);
            lowerLimit = structureStat(i).dvhStat.percDVH(lIx,:);
            upperLimit = structureStat(i).dvhStat.percDVH(hIx,:);
            confIn = structureStat(i).percentiles(hIx) - structureStat(i).percentiles(lIx);
            confName = ['C', num2str(round(confIn * 100,0))];
            h(hIx) = shadowPlot(doseGrid, lowerLimit, upperLimit, colors(j,:), confName, 1);
        end

        ylim([0 100]);
        if strcmp(pln.bioParam, 'RBExDose')
            xlabel('Dose RBE x [Gy]');
        else
            xlabel('Dose [Gy]');
        end
        ylabel('Volume [%]');
        lh = legend('show','Location','northeastoutside');
        uistack(h(2), 'top')
        uistack(h(1), 'top')
        labels = get(legend(), 'String');
        neworder = numel(labels):-1:1;
        plots = flipud(get(gca, 'children'));

        % Now re-create the legend
        legend(plots(neworder), labels(neworder))

        drawnow;
        hold off;
        cleanfigure();
        filename{i}.DVH = regexprep([cst{i,2},'_DVH.tex'], '\s+', '');
        matlab2tikz(fullfile(outputPath,'structures',filename{i}.DVH),'showInfo', false, 'width', '\figW', 'height', '\figH', 'extraAxisOptions', 'reverse legend');
        close

        % QI
        fprintf([num2str(i),'\n']);
        qiTable = fullNomQiTable(i,1:8);
        qiTable = [qiTable; structureStat(i).qiStat(:,1:8)];
        qiTable.Properties.RowNames{1} = 'nominal';
        input.data = qiTable;
        input.dataFormat = {'%.2f'};
        input.tableBorders = 0;
        input.booktabs = 1;
        input.tableEnvironment = 0;
        input.tableRowLabel = ['scen-', cst{i,2}];

        latex = latexTable(input);

        % save LaTex code as file
        filename{i}.QI = regexprep([cst{i,2},'_QI.tex'], '\s+', '');
        fid=fopen(fullfile(outputPath,'structures',filename{i}.QI),'w');
        [nrows,ncols] = size(latex);
        for row = 1:nrows
            fprintf(fid,'%s\n',latex{row,:});
        end
        fclose(fid);        
    end
end

% write them to structureWrapper
counter = 0;
line = cell(0);
for i = 1:size(cst,1)
    if cst{i,5}.Visible == true
        counter = counter + 1;
        if counter ~= 1
            line =  [line; '\newpage'];
        end
        line =  [line; ['\subsection{', cst{i,2}, '}']];
        line =  [line; '\begin{center}'];
        line =  [line; ['\input{', regexprep(fullfile(relativePath,filename{i}.DVH),'\','/'), '}']];
        line =  [line; '\end{center}'];
        line =  [line; ['\input{', regexprep(fullfile(relativePath,filename{i}.QI),'\','/'), '}']];
    end
end

fid = fopen(fullfile(outputPath,'structureWrapper.tex'),'w');
for i = 1:numel(line)
    text = regexprep(line{i},'\','\\\');
    fprintf(fid,text);
    fprintf(fid,'\n');
end
fclose(fid);


%% clean up
close all

end
