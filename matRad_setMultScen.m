function multScen = matRad_setMultScen(uIn)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_setPlanUncertainties function creates individual treatment planning scenarios
% 
% call
%   multScen = matRad_setMultScen(multScen)
%
% input
%   multScen:      multiple scenario struct
%
% output
%   multScen:      multiple scenario struct
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
  
if uIn.includeNomScen
   nomScen = 0;
else
   nomScen = [];
end

%% set isoCenter shifts
switch uIn.shiftGenType
    case 'equidistant'
        switch uIn.shiftGen1DIsotropy
            case '+-'
                % create grid vectors
                isoShiftVec{1} = [0 linspace(-uIn.shiftSize(1), uIn.shiftSize(1), uIn.numOfShiftScen(1))];
                isoShiftVec{2} = [0 linspace(-uIn.shiftSize(2), uIn.shiftSize(2), uIn.numOfShiftScen(2))];
                isoShiftVec{3} = [0 linspace(-uIn.shiftSize(3), uIn.shiftSize(3), uIn.numOfShiftScen(3))];
            case '+'
                isoShiftVec{1} = [0 linspace(0, uIn.shiftSize(1), uIn.numOfShiftScen(1))];
                isoShiftVec{2} = [0 linspace(0, uIn.shiftSize(2), uIn.numOfShiftScen(2))];
                isoShiftVec{3} = [0 linspace(0, uIn.shiftSize(3), uIn.numOfShiftScen(3))];        
            case '-'
                isoShiftVec{1} = [0 linspace(-uIn.shiftSize(1), 0, uIn.numOfShiftScen(1))];
                isoShiftVec{2} = [0 linspace(-uIn.shiftSize(2), 0, uIn.numOfShiftScen(2))];
                isoShiftVec{3} = [0 linspace(-uIn.shiftSize(3), 0, uIn.numOfShiftScen(3))];
        end
    case 'sampled'
        fprintf('sampled shifts only +- \n')
        % std is a 1 x 3 vector
        std = uIn.shiftSize;
        % mean (parameter)
        meanP = zeros(1,3);
        rng('shuffle');
        isoShiftVec{1} = [0 std(1) .* randn(1, uIn.numOfShiftScen(1)) + meanP(1)];
        rng('shuffle');
        isoShiftVec{2} = [0 std(2) .* randn(1, uIn.numOfShiftScen(2)) + meanP(2)];
        rng('shuffle');
        isoShiftVec{3} = [0 std(3) .* randn(1, uIn.numOfShiftScen(3)) + meanP(3)];     
    otherwise
        matRad_dispToConsole('did not expect that','error');
end

% create subMask for isoShifts as in the set masks part
% mask is a 'subClass' of the actual scenario mask
numIso(1) = numel(isoShiftVec{1});
numIso(2) = numel(isoShiftVec{2});
numIso(3) = numel(isoShiftVec{3});

scenMaskIso = false(numIso(1), numIso(2), numIso(3));

switch uIn.shiftCombType
    case 'individual'
        scenMaskIso(:,1,1) = true; % x shifts
        scenMaskIso(1,:,1) = true; % y shifts
        scenMaskIso(1,1,:) = true; % z shifts
    case 'permuted'
        scenMaskIso(:,:,:) = true;
    case 'combined'
        % determine that matrix is cubic
        if isequal(numIso(1), numIso(2), numIso(3))
            for i = 1:numIso(1)
                scenMaskIso(i,i,i) = true;
            end
        else
            uIn.shiftCombType = 'individual';
            matRad_dispToConsole('Numnber of isoShifts in every direction has to be equal in order to perform direct combination. Performing individually instead.',[],'warning');
            % call the function itself to get a working combination
            [multScen] = matRad_setUnc(uIn);
        end
    otherwise
        matRad_dispToConsole('Uncaught exception. Probably TYPO.','error');
end

% create list of increasing integers with referenced scenario
[xIso, yIso, zIso] = ind2sub(size(scenMaskIso),find(scenMaskIso));

matchMaskIso = cell(numel(xIso),2);
for i = 1:numel(xIso)
    matchMaskIso{i,1} = i;
    matchMaskIso{i,2} = [xIso(i) yIso(i) zIso(i)];
end

% create isoShift vectore based on the matrix and matching
isoShift = zeros(size(matchMaskIso,1),3);
if numel(isoShiftVec{1}) + numel(isoShiftVec{2}) + numel(isoShiftVec{3}) > 0
   
   for i = 1:size(matchMaskIso,1)
       matchPos = num2cell(matchMaskIso{i,2});
       
       if ~isequal([isoShiftVec{1}(matchPos{1}) isoShiftVec{2}(matchPos{2}) isoShiftVec{3}(matchPos{3})],[0 0 0]) || i == 1
             isoShift(i,:) = [isoShiftVec{1}(matchPos{1}) isoShiftVec{2}(matchPos{2}) isoShiftVec{3}(matchPos{3})] * ...
                             scenMaskIso(matchPos{:});
       end
   end
end

if isempty(isoShift)
   isoShift       = [0 0 0];
end

if ~uIn.includeNomScen
    isoShift = isoShift(2:end,:);
end
numOfShiftScen = size(isoShift,1);


%% range
switch uIn.rangeGenType
    case 'equidistant'
        relRangeShift = [nomScen linspace(-uIn.maxRelRangeShift, uIn.maxRelRangeShift, uIn.numOfRangeShiftScen)];
        absRangeShift = [nomScen linspace(-uIn.maxAbsRangeShift, uIn.maxAbsRangeShift, uIn.numOfRangeShiftScen)];
    case 'sampled'
        % relRange
        std = uIn.maxRelRangeShift;
        meanP = 0;
        rng('shuffle');
        relRangeShift = [nomScen std .* randn(1, uIn.numOfRangeShiftScen) + meanP];
        % absRange
        std = uIn.maxAbsRangeShift;
        meanP = 0;
        rng('shuffle');
        absRangeShift = [nomScen std .* randn(1, uIn.numOfRangeShiftScen) + meanP];
    otherwise
        matRad_dispToConsole('Not a valid type of generating data.','error');
end
numOfRelRangeShift = numel(relRangeShift);
numOfAbsRangeShift = numel(absRangeShift);

if ~isequal(numOfRelRangeShift, numOfAbsRangeShift)
    matRad_dispToConsole('Number of relative and absolute range shifts must not differ.','error');
else
    numOfRangeShiftScen = numOfRelRangeShift;
end

relRangeShift        = (relRangeShift./100);

% check if absolute and range error scenarios should be combined
switch uIn.rangeCombType    
   case 'individual'
      rangeShift = zeros(length(relRangeShift)*2 ,2);
      for i = 1:length(relRangeShift)
         rangeShift((i * 2)-1,1)  = absRangeShift(i);
         rangeShift((i * 2)  ,2)  = relRangeShift(i);
      end
      
      numOfRangeShiftScen = (numOfRelRangeShift*2)-1;
      absRangeShift       = [absRangeShift zeros(1,length(relRangeShift(relRangeShift~=0)))];
      relRangeShift       = [zeros(1,length(relRangeShift)) relRangeShift(relRangeShift~=0)] ;
      
   case 'combined' 
      rangeShift = zeros(length(relRangeShift),2);
      for i = 1:1:length(relRangeShift) 
         rangeShift(i,1) = absRangeShift(i);
         rangeShift(i,2) = relRangeShift(i);
      end
      
   otherwise
      
end

% check if the existence shift scenarios
 if sum(sum(~any(isoShift))) > 0 && ~uIn.includeNomScen
    isoShift = [];
 end

% combine setup and range scenarios according to scenCombType
switch uIn.scenCombType
      
    case 'individual' % combine setup and range scenarios individually

       % range errors should come first
       if uIn.includeNomScen
          scenForProb                                = zeros(size(isoShift,1)-1 + size(rangeShift,1),5);
          scenForProb(1:size(rangeShift,1),4:5)      = rangeShift;
          scenForProb(size(rangeShift,1)+1:end,1:3)  = isoShift(2:end,:);
       else
           scenForProb                               = zeros(size(isoShift,1) + size(rangeShift,1),5);
           scenForProb(1:size(rangeShift,1),4:5)     = rangeShift;
           scenForProb(size(rangeShift,1)+1:end,1:3) = isoShift;
       end
       
    case 'permuted'

        scenForProb  = zeros(size(isoShift,1) * size(rangeShift,1),5);
        Cnt = 1;
        for i = 1:size(isoShift,1)
           for j = 1:size(rangeShift,1)
              scenForProb(Cnt,:)   = [isoShift(i,:) rangeShift(j,:)];
              Cnt = Cnt + 1;
           end
        end

    case 'combined'

       if size(isoShift,1) == size(rangeShift,1)  && numOfShiftScen > 0 && numOfRangeShiftScen > 0
        scenForProb            = zeros(size(isoShift,1),5);
        scenForProb(1:end,1:3) = isoShift;
        scenForProb(1:end,4:5) = rangeShift; 

       else
           matRad_dispToConsole('number of setup and range scenarios MUST be the same \n',[],'warning');
           uIn.scenCombType = 'individual';
           multScen         = matRad_setMultScen(uIn);
           scenForProb      = multScen.scenForProb;
       end
end

% if no shift scenarios are defined isoShift needs to be set to [0 0 0] to ensure a proper dose calc
if isempty(isoShift)
   isoShift = [0 0 0];
end

% sanity check
UniqueRowScenForProb = unique(scenForProb,'rows');

if size(UniqueRowScenForProb,1) ~= size(scenForProb,1) && size(UniqueRowScenForProb,1)>1
     matRad_dispToConsole('Some scenarios seem to be defined multiple times',[],'warning');
end

multScen.scenForProb = scenForProb;

%% set masks

% 1st dim = rows: isoShift, 2nd dim: relRangeShift, 3rd dim: absRangeShift
scenMask        = false(uIn.numOfCtScen, numOfShiftScen, numOfRangeShiftScen);
scenMask(:,1,1) = true; % ct scenarios

% switch between combination modes here
% only makes scence when numOfShiftScen>0 and numOfRangeShiftScen>0;
if numOfShiftScen > 0 && numOfRangeShiftScen > 0
   switch uIn.scenCombType
       case 'individual'
          
            % get all setup scenarios
           [~,ixUnq] = unique(scenForProb(:,1:3),'rows','stable');
          
           scenMask  = false(uIn.numOfCtScen, length(ixUnq), numOfRangeShiftScen);
          
           scenMask(1,:,1) = true; % iso shift scenarios
           scenMask(1,1,:) = true; % range shift scenarios
           
       case 'permuted'
          
           scenMask(:,:,:) = true;
           
       case 'combined'
          
           % check if scenForProb is cubic (ignore ct scen) - if so fill diagonal
           if isequal(numOfShiftScen, numOfRangeShiftScen)
               for i = 1:size(scenForProb, 1)
                   scenMask(1,i,i) = true;
               end
           else
               uIn.shiftCombType = 'individual';
               matRad_dispToConsole('number of setup and range scenarios MUST be the same \n',[],'warning');
               scenMask = multScen.scenMask;
           end
   end
end

% create linearalized mask where the i row points to the indexes of scenMask
[x{1}, x{2}, x{3}]           = ind2sub(size(scenMask),find(scenMask));
linearMask                   = cell2mat(x);

% get number of scenarios
totalNumScen                 = size(scenForProb, 1);

% write to scenario information variable
multScen.numOfCtScen         = uIn.numOfCtScen;

multScen.isoShift            = isoShift;
multScen.shiftCombType       = uIn.shiftCombType;
multScen.shiftGen1DIsotropy  = uIn.shiftGen1DIsotropy;
multScen.numOfShiftScen      = numOfShiftScen;

multScen.relRangeShift       = relRangeShift;
multScen.absRangeShift       = absRangeShift;
multScen.numOfRangeShiftScen = numOfRangeShiftScen;
multScen.rangeCombType       = uIn.rangeCombType;

multScen.scenMask            = scenMask;
multScen.linearMask          = linearMask;
multScen.numOfScen           = totalNumScen;

multScen.rangeRelSD          = uIn.rangeRelSD;
multScen.rangeAbsSD          = uIn.rangeAbsSD;
multScen.shiftSD             = uIn.shiftSD;  

multScen.scenCombType        = uIn.scenCombType;

% save original user input
multScen.userInput = uIn;



