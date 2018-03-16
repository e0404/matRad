function accDose = calc4dDose_loop(ct, pln, dij, stf, cst, resultGUI, FileName, GTVName, OARName, count, MOTION, phaseTimeDev)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculation of accumulated dose for large number of Lmdout files  
%
% please check: dose accumulation DDM or EMT and correct VF?!?
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 12
    MOTION = 'linear';
    phaseTimeDev = 0;
end

addpath('D:\Matrad\')
addpath('D:\Matrad\4Ddose')

resultGUI.bioParam = pln.bioParam;

c=1;
Ndoses = 1;
maxnumberoffset = 5;
maxnumbermotionperiod = 5;
maxnumbermotionvariations = 1;
%offsetA = [0 1.5 3 2 8 6 4 1 7 3 9 4.5 6.5 2.5 8.5 0.5 7.5 1.5 3.5 0.2 7 1.8 4.5 0.4 10 0.3 4.1 2.1 9.4 6.9 5.5];
%motionperiodA = [5 3 8 10 4 6 5 7 5.5 6.5 3.5 9.5 4.5 7.5 9.5];

NumOfPhases = size(ct.cube,2);
if isfield(ct, 'dvf')
if size(ct.dvf,2) < size(ct.cube,2)
    NumOfPhases = size(ct.dvf,2);
end
end

while c<=count
    disp(['Delivery file ', num2str(c), ' von ', num2str(count)]);
    if(count ==1 )
       FileName_lmdout = FileName; 
    else
        n = num2str(c);
        FileName_lmdout = [FileName '_' n];
    end
    
    % reads in PB XML Plan and result of dose delivery simulation and creates
    % delievery struct    
    delivery = matRad_readLmdout(dij, stf, FileName, FileName_lmdout);
    delivery(1).NumOfPhases = NumOfPhases;
    
    % dose in each CT phase is calculated
    for offset = 1:maxnumberoffset
        disp(['offset ', num2str(offset), ' von ', num2str(maxnumberoffset)]);
        for motionperiod = 1:maxnumbermotionperiod
            disp(['period ', num2str(motionperiod), ' von ', num2str(maxnumbermotionperiod)]);
            for motionvariation = 1: maxnumbermotionvariations
                 disp(['variation ', num2str(motionvariation), ' von ', num2str(maxnumbermotionvariations)]);
                 
                 % Breathing period is assumed as a normal distribution
                 % centered at 6 sec
                 delivery(1).motionperiod = 6+ randn(1); % breathing period  3-10 ? 5? %motionperiodA(motionperiod);                 
                 delivery(1).offset = NumOfPhases*rand(1); %offsetA(offset); 
                 
    %[resultGUI, delivery] = matRad_calcPhaseDose(resultGUI, dij, delivery);  
    [resultGUI, delivery] = matRad_calcPhaseDoseMatrix(resultGUI, dij, delivery, MOTION, phaseTimeDev);  %Matrix
    
    % dose accumulation
    resultGUI = matRad_doseAcc(ct, resultGUI, cst, 'DDM');  %acc Methods: 'EMT' 'DDM'
    
    if  strcmp(pln.bioOptimization,'none')
        accDose.D{Ndoses} = resultGUI.accDose;
    else
        accDose.D{Ndoses} = resultGUI.accRBExDose;    
    end
    
    accDose.offset{Ndoses} = delivery(1).offset;
    accDose.motionperiod{Ndoses} = delivery(1).motionperiod;
    
    Ndoses = Ndoses+1;
    
        end
        end
    end
    c=c+1;
end


%plot
vIsoCenter = round(pln.isoCenter./[ct.resolution.x ct.resolution.y ct.resolution.z]);
slice = vIsoCenter(3);

fignr = size(accDose.D);
fignr = fignr(2);

c= ceil(fignr/3);
r = ceil(fignr/c);



%%
%DVH Statistiken
%für GTV
for i=1:size(cst,1)
        if strcmp(cst{i,2},GTVName)
            indices     = cst{i,4}{1}; % GTV %cst{5,4}{1};  %GTV für Liver007  
        end
end

count = count*maxnumberoffset*maxnumbermotionperiod*maxnumbermotionvariations;
refVol = [5 95 1 99]; %5% und 95%
numOfVoxels = numel(indices);
minD99 = 10000;
maxD1 = 0;
for f=1:count
    
    % get Dose, dose is sorted to simplify calculations
    relevantDose = accDose.D{f};
    doseInVoi    = sort(accDose.D{f}(indices));
    
    %refGy = round([0.4 0.6 0.8] * max(relevantDose(:)) * 10)/10;
    refGy = round(linspace(1,round(max(relevantDose(:))),3));
    
    DX = @(x) doseInVoi(ceil((100-x)*0.01*numOfVoxels));
    VX = @(x) numel(doseInVoi(doseInVoi >= x)) / numOfVoxels;
      
    for runDX = 1:numel(refVol)
       QI.(strcat('D',num2str(refVol(runDX)))) = DX(refVol(runDX));          
    end    
    
    accDose.D95(f) = QI.D95;
        
    accDose.HI(f) = (DX(2)/DX(98));
    
    if(QI.D99 < minD99)
        minD99 = QI.D99;
    end
    if(QI.D1 > maxD1)
        maxD1 = QI.D1;
    end        
end

accDose.meanD95 = mean(accDose.D95);
accDose.stdD95 = std(accDose.D95);

accDose.meanHI = mean(accDose.HI);
accDose.stdHI = std(accDose.HI);

accDose.minD99 = minD99;
accDose.maxD1 = maxD1;


%statisch
f=f+1;
 if  strcmp(pln.bioOptimization,'none')
relevantDose = resultGUI.physicalDose;
doseInVoi    = sort(resultGUI.physicalDose(indices));
 else
   relevantDose = resultGUI.RBExDose;
 doseInVoi    = sort(resultGUI.RBExDose(indices));
 end
 refGy = round(linspace(1,round(max(relevantDose(:))),3));
    
    DX = @(x) doseInVoi(ceil((100-x)*0.01*numOfVoxels));
    VX = @(x) numel(doseInVoi(doseInVoi >= x)) / numOfVoxels;
  
     
    accDose.D95(f) = DX(95); 
    accDose.HI(f) = (DX(2)/DX(98));
      
%%
     
     

 %1 DVH
figure
n=1000;
 if  strcmp(pln.bioOptimization,'none')
dvhPointsOpt = linspace(0,max(resultGUI.finalDose(:))*1.05,n);
dvh       = NaN * ones(1,n);
%indices     = cst{5,4}{1};  %GTV
numOfVoxels = numel(indices);
doseInVoi   = resultGUI.finalDose(indices);
 else
    dvhPointsOpt = linspace(0,max(resultGUI.RBExDose(:))*1.05,n);
dvh       = NaN * ones(1,n);
%indices     = cst{5,4}{1};  %GTV
numOfVoxels = numel(indices);
doseInVoi   = resultGUI.RBExDose(indices); 
 end
for j = 1:n
    dvh(j) = sum(doseInVoi > dvhPointsOpt(j));
end
dvhOpt = dvh ./ numOfVoxels * 100;
  
    %plot(dvhPointsOpt, dvhOpt, 'LineWidth', 2); hold on 
for f=1:count
    dvhPoints = linspace(0,max(accDose.D{f}(:))*1.05,n);
    dvh       = NaN * ones(1,n);

   % indices     = cst{5,4}{1};  %GTV
    numOfVoxels = numel(indices);
    doseInVoi   = accDose.D{f}(indices);

    for j = 1:n
            dvh(j) = sum(doseInVoi > dvhPoints(j));
    end

    dvh = dvh ./ numOfVoxels * 100;
    plot(dvhPoints, dvh, 'LineWidth', 1.5);hold on
end
plot(dvhPointsOpt, dvhOpt, 'LineWidth', 2, 'Color', 'r');
%title(['GTV: max D5% = ' num2str(maxD5) '  min D95% = ' num2str(minD95)])       
title(['DVH ' GTVName])
xlabel('RBExDose [Gy(RBE)]','FontSize',14)
ylabel('Volume [%]', 'FontSize',14)
set(gca,'FontSize',14) 


%Liver

for i=1:size(cst,1)
        if strcmp(cst{i,2},OARName)

indices     = cst{i,4}{1};

        end
end

refVol = [1 50]; 
numOfVoxels = numel(indices);

for f=1:count
    
    % get Dose, dose is sorted to simplify calculations
    relevantDose = accDose.D{f};
    doseInVoi    = sort(accDose.D{f}(indices));
    
    %refGy = round([0.4 0.6 0.8] * max(relevantDose(:)) * 10)/10;
    refGy = round(linspace(1,round(max(relevantDose(:))),3));
    
    DX = @(x) doseInVoi(ceil((100-x)*0.01*numOfVoxels));
    VX = @(x) numel(doseInVoi(doseInVoi >= x)) / numOfVoxels;
      
    for runDX = 1:numel(refVol)
       QI.(strcat('D',num2str(refVol(runDX)))) = DX(refVol(runDX));          
    end    
    
    accDose.D1Liver(f) = QI.D1;
        
    accDose.D50Liver(f) = QI.D50;
    accDose.meanLiver(f) =  mean(doseInVoi);    
end

accDose.meanD1Liver = mean(accDose.D1Liver);
accDose.stdDD1Liver = std(accDose.D1Liver);

accDose.meanmeanLiver = mean(accDose.meanLiver);
accDose.stdmeanLiver = std(accDose.meanLiver);

%statisch 
f=f+1;
 if  strcmp(pln.bioOptimization,'none')
relevantDose = resultGUI.physicalDose;
doseInVoi    = sort(resultGUI.physicalDose(indices));
 else
   relevantDose = resultGUI.RBExDose;
doseInVoi    = sort(resultGUI.RBExDose(indices));
 end
 refGy = round(linspace(1,round(max(relevantDose(:))),3));
    
    DX = @(x) doseInVoi(ceil((100-x)*0.01*numOfVoxels));
    VX = @(x) numel(doseInVoi(doseInVoi >= x)) / numOfVoxels;
  
     
    accDose.D1Liver(f) = QI.D1;
    accDose.meanLiver(f) = mean(doseInVoi);


%1 DVH OAR
figure
n=1000;
 if  strcmp(pln.bioOptimization,'none')
dvhPointsOpt = linspace(0,max(resultGUI.finalDose(:))*1.05,n);
dvh       = NaN * ones(1,n);
numOfVoxels = numel(indices);
doseInVoi   = resultGUI.finalDose(indices);
 else
    dvhPointsOpt = linspace(0,max(resultGUI.RBExDose(:))*1.05,n);
dvh       = NaN * ones(1,n);
numOfVoxels = numel(indices);
doseInVoi   = resultGUI.RBExDose(indices); 
 end
for j = 1:n
    dvh(j) = sum(doseInVoi > dvhPointsOpt(j));
end
dvhOpt = dvh ./ numOfVoxels * 100;
  
for f=1:count
    dvhPoints = linspace(0,max(accDose.D{f}(:))*1.05,n);
    dvh       = NaN * ones(1,n);

    numOfVoxels = numel(indices);
    doseInVoi   = accDose.D{f}(indices);

    for j = 1:n
            dvh(j) = sum(doseInVoi > dvhPoints(j));
    end

    dvh = dvh ./ numOfVoxels * 100;
    plot(dvhPoints, dvh, 'LineWidth', 1.5);hold on
end
plot(dvhPointsOpt, dvhOpt, 'LineWidth', 2, 'Color', 'r');   
title(['DVH ' OARName])
xlabel('RBExDose [Gy(RBE)]','FontSize',14)
ylabel('Volume [%]', 'FontSize',14)
set(gca,'FontSize',14) 

