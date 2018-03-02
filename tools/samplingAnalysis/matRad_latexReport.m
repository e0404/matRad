function param = matRad_latexReport(ct, cst, pln, nominalScenario, structureStat, doseStat, sampDose, listOfQI, param)
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

function text = parseFromDicom(dicomStruct, field, default)
    if ~exist('default', 'var') || isempty(default)
        default = 'N.A.';
    end
    if isfield(dicomStruct, field)
        text = dicomStruct.(field);
    else
        text = default;
    end    
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

%% check whether all qi are available
if  ~all(ismember(listOfQI, fieldnames(nominalScenario.qi)))
    error('Not all QI are available. Possible TYPO in at least one of them.');
end
[~, ~, ixOfQIinNominal] = intersect(listOfQI, fieldnames(nominalScenario.qi), 'stable');

% check param
optionalFields = {'computationTime'};
for i=1:numel(optionalFields)
    if ~isfield(param, optionalFields{i})
        param.(optionalFields{i}) = NaN;
    end
end
        
%% insert standard patient information
% issue warning for insufficient number of scenarios
if exist('param','var') && isfield(param,'sufficientStatistics') && ~param.sufficientStatistics
    warnMessage = 'Insufficient statistics. Handle with care.';
else
    warnMessage = '';
end

if isfield(ct, 'dicomInfo')
    if isfield(ct.dicomInfo, 'PatientName')
        patientInformation.firstName = parseFromDicom(ct.dicomInfo.PatientName, 'GivenName');
        patientInformation.lastName = parseFromDicom(ct.dicomInfo.PatientName, 'FamilyName');
    end
    %%% in future release, patient sex and patient id might be moved in dicomInfo
    patientInformation.sex = parseFromDicom(ct.dicomInfo, 'PatientSex');
    patientInformation.patientID = parseFromDicom(ct.dicomInfo, 'PatientID');
else
    patientInformation.firstName = 'N.A.';
    patientInformation.lastName = 'N.A.';
    patientInformation.sex = 'N.A.';
    patientInformation.patientID = 'N.A.';
end

%%% when id and sex moved to dicomInfo this get's obsolete
try
    patientInformation.sex = ct.dicomMeta.PatientSex;
catch
end
try
    patientInformation.patientID = ct.dicomMeta.PatientID;
catch
end

% import plan information
planInformation.gantryAngles = num2str(pln.gantryAngles);
planInformation.couchAngles = num2str(pln.couchAngles);
planInformation.modality = pln.radiationMode;

line = cell(0);
line =  [line; '\newcommand{\warnMessage}{',warnMessage,'}'];
line =  [line; '\newcommand{\patientFirstName}{',patientInformation.firstName,'}'];
line =  [line; '\newcommand{\patientLastName}{',patientInformation.lastName,'}'];
line =  [line; '\newcommand{\patientSex}{',patientInformation.sex,'}'];
line =  [line; '\newcommand{\patientID}{',patientInformation.patientID,'}'];
line =  [line; '\newcommand{\operator}{',param.operator,'}'];

line =  [line; '\newcommand{\reportGenerationDate}{\today}'];
line =  [line; '\newcommand{\computationTime}{', num2str(round(param.computationTime / 3600,2)), ' h}'];

line =  [line; '\newcommand{\planGantryAngles}{',planInformation.gantryAngles,'}'];
line =  [line; '\newcommand{\planCouchAngles}{',planInformation.couchAngles,'}'];
line =  [line; '\newcommand{\planRadiationModality}{',planInformation.modality,'}'];

fid = fopen(fullfile(outputPath,'patientInformation.tex'),'w');
for i = 1:numel(line)
    text = regexprep(line{i},{'\','_'},{'\\\','-'});
    fprintf(fid,text);
    fprintf(fid,'\n');
end
fclose(fid);

%% analysis parameters
line = cell(0);

% raw input
multScen = pln.multScen;
% shift parameters
line =  [line; '\newcommand{\numOfShiftScen}{', num2str(multScen.numOfShiftScen), '}'];
line =  [line; '\newcommand{\shiftSize}{', num2str(max(multScen.isoShift)), '}'];
line =  [line; '\newcommand{\shiftGenType}{', num2str(multScen.shiftGenType), '}'];
line =  [line; '\newcommand{\shiftCombType}{', num2str(multScen.shiftCombType), '}'];

% range parameters
line =  [line; '\newcommand{\numOfRangeShiftScen}{', num2str(multScen.numOfRangeShiftScen), '}'];
line =  [line; '\newcommand{\maxAbsRangeShift}{', num2str(max(multScen.absRangeShift)), '}'];
line =  [line; '\newcommand{\maxRelRangeShift}{', num2str(max(multScen.relRangeShift)), '}'];
line =  [line; '\newcommand{\rangeCombType}{', num2str(multScen.rangeCombType), '}'];
line =  [line; '\newcommand{\rangeGenType}{', num2str(multScen.rangeGenType), '}'];
line =  [line; '\newcommand{\scenCombType}{', num2str(multScen.scenCombType), '}'];

% gamma analysis parameters
line =  [line; '\newcommand{\gammaDoseAgreement}{', num2str(doseStat.gammaAnalysis.doseAgreement), '}'];
line =  [line; '\newcommand{\gammaDistAgreement}{', num2str(doseStat.gammaAnalysis.distAgreement), '}'];

if pln.multScen.numOfCtScen <= 1
    line =  [line; '\newcommand{\ctScen}{false}'];
else
    line =  [line; '\newcommand{\ctScen}{true}'];
end

if pln.multScen.numOfRangeShiftScen <= 1
    line =  [line; '\newcommand{\rangeScen}{false}'];
else
    line =  [line; '\newcommand{\rangeScen}{true}'];
end

if pln.multScen.numOfShiftScen <= 1
    line =  [line; '\newcommand{\shiftScen}{false}'];
else
    line =  [line; '\newcommand{\shiftScen}{true}'];
end

line =  [line; '\newcommand{\shiftSD}{', num2str(pln.multScen.shiftSD), '}'];
line =  [line; '\newcommand{\rangeAbsSD}{', num2str(pln.multScen.rangeAbsSD), '}'];
line =  [line; '\newcommand{\rangeRelSD}{', num2str(pln.multScen.rangeRelSD), '}'];

fid = fopen(fullfile(outputPath,'uncertaintyParameters.tex'),'w');
for i = 1:numel(line)
    text = regexprep(line{i},'\','\\\');
    fprintf(fid,text);
    fprintf(fid,'\n');
end
fclose(fid);


%% plot isocentre slices (nominal, mean, std)
doseCube = nominalScenario.(pln.bioParam.quantityVis);

doseWindow  = [0 1.1 * max(doseCube(:))];
stdWindow   = [0 1.1 * max(doseStat.stdCubeW(:))];
gammaWindow = [0 1.1 * max(doseStat.gammaAnalysis.gammaCube(:))];

for plane=1:3
    switch plane 
        case 1
            slice = round(pln.isoCenter(1,2) / ct.resolution.x,0);
        case 2
            slice = round(pln.isoCenter(1,1) / ct.resolution.y,0);
        case 3
            slice = round(pln.isoCenter(1,3) / ct.resolution.z,0);
    end
    colors = colorcube(size(cst,1));
    for cubesToPlot = 1:3
        figure; ax = gca;
        
        doseCube = nominalScenario.(pln.bioParam.quantityVis); 
        
        if cubesToPlot == 1 
            
            if isfield(nominalScenario,'RBExD')              
              colorMapLabel = 'RBExDose [Gy(RBE)]';
            else
                colorMapLabel = 'physical Dose [Gy]';
            end
            fileSuffix = 'nominal';
            matRad_plotSliceWrapper(ax,ct,cst,1,doseCube,plane,slice,[],[],colors,[],[],[],[],colorMapLabel);
            
        elseif cubesToPlot == 2
            
            doseCube = doseStat.gammaAnalysis.gammaCube;
            colorMapLabel = 'gamma index';
            fileSuffix = 'gamma';
            gammaColormap = matRad_getColormap('gammaIndex');
            matRad_plotSliceWrapper(ax,ct,cst,1,doseCube,plane,slice,[],[],colors,gammaColormap,[0 2],[],[],colorMapLabel);
            
        elseif cubesToPlot == 3
            
            if isfield(nominalScenario,'RBExD')
                colorMapLabel = 'Standard deviation [Gy(RBE)]';
            else
                colorMapLabel = 'Standard deviation [Gy]';
            end
            doseCube = doseStat.stdCubeW;
            fileSuffix = 'stdW';
            matRad_plotSliceWrapper(ax,ct,cst,1,doseCube,plane,slice,[],[],colors,[],[],[],[],colorMapLabel);
            
        end
        drawnow();
        cleanfigure();          
        matlab2tikz(fullfile(outputPath,['isoSlicePlane', num2str(plane), '_', fileSuffix, '.tex']), 'relativeDataPath', 'data', 'showInfo', false, 'width', '\figW');
        close
    end
end

%% gaussian orbit sample animation
if exist('matRad_getGaussianOrbitSamples','file') == 2
    param.confidenceValue = 0.5;

    slice = round(pln.isoCenter(1,plane) / ct.resolution.z,0);            
    outPath = fullfile(outputPath, 'frames');
    if isfield(nominalScenario,'RBExD')
        legendColorbar = 'RBExDose [Gy(RBE)]';
    else
        legendColorbar = 'physical Dose [Gy]';
    end
    % any(w == 0) is not allowed, due to numerical reasons use insignificant w, for weights which are numerically zero
    w = pln.multScen.scenProb;
    w(w == 0) = eps(class(w));
    matRad_createAnimationForLatexReport(param.confidenceValue, ct, cst, slice, doseStat.meanCubeW, sampDose, w, pln.subIx, outPath, legendColorbar);

    line = cell(0);
    line =  [line; '\newcommand{\framerate}{24}'];
    line =  [line; '\newcommand{\firstframe}{1}'];
    line =  [line; '\newcommand{\lastframe}{120}'];

    fid = fopen(fullfile(outputPath,'parameters.tex'),'w');
    for i = 1:numel(line)
        text = regexprep(line{i},{'\','_'},{'\\\','-'});
        fprintf(fid,text);
        fprintf(fid,'\n');
    end
    fclose(fid);
end

%% plot nominal dvh and write qi table
% plot dvh
colors = jet(size(cst,1));
hold off;
for i = 1:size(cst,1)
    if cst{i,5}.Visible == true
        [y, argmin] = cutAtArgmin(nominalScenario.dvh(i).volumePoints);
        x = nominalScenario.dvh(i).doseGrid(1:argmin);        
        h(1) = plot(x,y,'LineWidth',2, 'Color', colors(i,:), 'DisplayName', cst{i,2});      
        ylim([0 100]);
        if strncmp(pln.bioParam.quantityVis,'RBExD',5)
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
fieldsInQi = fieldnames(nominalScenario.qi);
for i= 1:numel(fieldsInQi)
    nanRow.(fieldsInQi{i}) = NaN;
end

nomQi = nominalScenario.qi(1);
structName = {};
c = 1;
for i = 1:size(cst,1)
    if cst{i,5}.Visible == true && ~isempty(cst{i,4}{1})
        nomQi(c) = nominalScenario.qi(i);
        fullQi(i) = nominalScenario.qi(i);
        structName = [structName, cst{i,2}];
        c = c + 1;
    else
        fullQi(i) = nanRow;
    end
end



nomQiTable = struct2table(nomQi);
fullNomQiTable = struct2table(fullQi);
fullNomQiTable = fullNomQiTable(:,ixOfQIinNominal);
nomQiTable.Properties.RowNames = structName;
nomQiTable.Properties.VariableNames = regexprep(nomQiTable.Properties.VariableNames, '_', '');
input.data = nomQiTable(:,ixOfQIinNominal);
input.dataFormat = {'%.2f'};
input.booktabs = 1;
input.tableBorders = 0;
input.tableEnvironment = 0;

latex = latexTable(input);

% save LaTex code as file
filename{i}.QI = regexprep([cst{i,2},'_QI.tex'], '\s+', '');
fid=fopen(fullfile(outputPath,'nominalQI.tex'),'w');
[nrows, ~] = size(latex);
for row = 1:nrows
    fprintf(fid,'%s\n',latex{row,:});
end
fclose(fid);
clear latex
clear input

%% add DVH and QI
clear h
clear filename
% relative file path (relative to main.tex)
relativePath = fullfile('data','structures');

if strcmp(pln.bioParam.quantityVis, 'RBExD')
    labelDoseDVH = 'Dose RBE x [Gy]';
else
   labelDoseDVH = 'Dose [Gy]';
end
for i = 1:size(cst,1)
    if cst{i,5}.Visible == true
        [~, ~, ixOfQIinStruct] = intersect(listOfQI, fieldnames(structureStat(i).qiStat), 'stable');
        
        matRad_plotDVHBand(nominalScenario.dvh(i), structureStat(i), labelDoseDVH);
        cleanfigure();
        filename{i}.DVH = regexprep([cst{i,2},'_DVH.tex'], '\s+', '');
        matlab2tikz(fullfile(outputPath,'structures',filename{i}.DVH),'showInfo', false, 'width', '\figW', 'height', '\figH', 'extraAxisOptions', 'reverse legend');
        close
        
        % QI
        fprintf([num2str(i),'\n']);
        qiTable = fullNomQiTable(i,:);
        qiTable = [qiTable; structureStat(i).qiStat(:,ixOfQIinStruct)];
        qiTable.Properties.RowNames{1} = 'nominal';
        qiTable.Properties.VariableNames = regexprep(qiTable.Properties.VariableNames, '_', '');
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
        [nrows,~] = size(latex);
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
