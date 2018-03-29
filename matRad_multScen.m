classdef matRad_multScen 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  matRad_multScen 
%  This class creates all required biological model parameters according to
% a given radiation modatlity and a given bio model identifier string. 
%
% constructor call
%   pln.multScen = matRad_multScen(ct,TYPE);
%
%   e.g. pln.multScen = matRad_multScen(ct,'nomScen');
%
% input
%   ct:                 ct cube
%   TYPE:               string to denote a scenario creation method
%                       'nomScen'   create only the nominal scenario
%                       'wcScen'    create worst case scenarios
%                       'impScen'   create important/grid scenarios
%                       'rndScen'   create random scenarios
%
% output
%   bioParam:           matRad's bioParam structure containing information
%                       about the choosen biological model
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
    
    %% properties
    
    % public properties which can be changed outside this class
    properties               
             
        TYPE;                   % denotes the sampling type which cen be one of the following 
                                % 'nomScen'   create only the nominal scenario
                                % 'wcScen'    create worst case scenarios
                                % 'impScen'   create important/grid scenarios
                                % 'rndScen'   create random scenarios
        % ct scenarios
        numOfCtScen;            % number of imported ct scenarios                        
        
        % shift scenarios
        numOfShiftScen;         % number of shifts in x y and z direction       
        shiftSize;              % 3x1 vector to define maximal shift in [mm]  % (e.g. abdominal cases 5mm otherwise 3mm)
        shiftGenType;           % 'equidistant': equidistant shifts, 'sampled': sample shifts from normal distribution
        shiftCombType;          % 'individual':  no combination of shift scenarios;       number of shift scenarios is sum(multScen.numOfShiftScen) 
                                % 'combined':    combine shift scenarios;                 number of shift scenarios is numOfShiftScen
                                % 'permuted':    create every possible shift combination; number of shift scenarios is 8,27,64 .
                        
        % b) define range error scenarios                                                
        numOfRangeShiftScen;    % number of absolute and/or relative range scnearios. 
                                % if absolute and relative range scenarios are defined then rangeCombType defines the resulting number of range scenarios
        maxAbsRangeShift;       % maximum absolute over and undershoot in mm   
        maxRelRangeShift;       % maximum relative over and undershoot in % 
        rangeCombType;          % 'individual':  no combination of absolute and relative range scenarios
                                % 'combined':    combine absolute and relative range scenarios
        rangeGenType;           % 'equidistant': equidistant range shifts, 'sampled': sample range shifts from normal distribution
        scenCombType;           % 'individual':  no combination of range and setup scenarios, 
                                % 'combined':    combine range and setup scenarios if their scenario number is consistent 
                                % 'permuted':    create every possible combination of range and setup scenarios
                                                 
        includeNomScen;         % boolean to determine if the nominal scenario should be included
         
        % these parameters will be filled according to the choosen scenario type
        isoShift;
        relRangeShift;
        absRangeShift;
            
        totNumShiftScen;        % total number of shift scenarios in x,y and z direction
        totNumRangeScen;        % total number of range and absolute range scenarios
        totNumScen;             % total number of samples
        
        scenForProb;            % matrix for probability calculation - each row denotes one scenario  
        scenProb;               % probability of each scenario stored in a vector
        scenMask;
        linearMask;
        
        rangeRelSD  = 3.5;                % given in %
        rangeAbsSD  = 1;                  % given in [mm]   
        shiftSD     = [2.25 2.25 2.25];   % given in [mm]
        
    end
    
    % private properties which can only be changed inside matRad_multScen
    properties(SetAccess = private)
        
        DEFAULT_TYPE      = 'nomScen'; 
        DEFAULT_probDist  = 'normDist'; % 'normDist': normal probability distrubtion or 'equalProb' for uniform probability distribution  
        DEFAULT_TypeProp  = 'probBins'; % 'probBins': cumulative prob in bin or 'pointwise' for point probability  
        
    end
    
    % constant properties which are visible outside of matRad_multScen
    properties(Constant = true)
        
        AvailableScenCreationTYPE = {'nomScen','wcScen','impScen','rndScen'};
    end
    
    % constant private properties which are only visible within matRad_multScen
    properties(Constant = true, Access = private)
        
        % default parameter for each scenario creation TYPE
        
        % 'nomScen' default parameters for nominal scenario
        numOfShiftScen_nomScen            = [0 0 0];                       % number of shifts in x y and z direction
        shiftSize_nomScen                 = [0 0 0];                       % given in [mm]
        shiftGenType_nomScen              = 'equidistant';                 % equidistant: equidistant shifts,
        shiftCombType_nomScen             = 'individual';                  % individual:  no combination of shift scenarios;
        numOfRangeShiftScen_nomScen       = 0                              % number of absolute and/or relative range scnearios.
        maxAbsRangeShift_nomScen          = 0;                             % maximum absolute over and undershoot in mm 
        maxRelRangeShift_nomScen          = 0;                             % maximum relative over and undershoot in % 
        rangeCombType_nomScen             = 'combined';                    % combine absolute and relative range scenarios
        rangeGenType_nomScen              = 'equidistant';                 % equidistant: equidistant range shifts,
        scenCombType_nomScen              = 'individual';                  % individual:  no combination of range and setup scenarios, 
        includeNomScen_nomScen            = true;
                             
        % 'wcScen'  default parameters for  worst case scenarios
        numOfShiftScen_wcScen             = [2 2 2];                       % number of shifts in x y and z direction
        shiftSize_wcScen                  = [3 3 3];                       % given in [mm]
        shiftGenType_wcScen               = 'equidistant';                 % equidistant: equidistant shifts
        shiftCombType_wcScen              = 'individual';                  % individual:  no combination of shift scenarios
        numOfRangeShiftScen_wcScen        = 2;                             % number of absolute and/or relative range scnearios.
        maxAbsRangeShift_wcScen           = 1;                             % maximum absolute over and undershoot in mm 
        maxRelRangeShift_wcScen           = 3.5;                           % maximum relative over and undershoot in % 
        rangeCombType_wcScen              = 'combined';                    % combine absolute and relative range scenarios
        rangeGenType_wcScen               = 'equidistant';                 % equidistant: equidistant range shifts
        scenCombType_wcScen               = 'individual';                  % individual:  no combination of range and setup scenarios
        includeNomScen_wcScen             = true;                          % include nominal scenario
        
        % 'impScen'   create important/grid scenarios                       
        numOfShiftScen_impScen            = [0 0 0];                       % number of shifts in x y and z direction
        shiftSize_impScen                 = [0 0 0];                       % given in [mm]
        shiftGenType_impScen              = 'equidistant';                 % equidistant: equidistant shifts
        shiftCombType_impScen             = 'individual';                  % individual:  no combination of shift scenarios
        numOfRangeShiftScen_impScen       = 20;                            % number of absolute and/or relative range scnearios
        maxAbsRangeShift_impScen          = 1;                             % maximum absolute over and undershoot in mm 
        maxRelRangeShift_impScen          = 3.5;                           % maximum relative over and undershoot in % 
        rangeCombType_impScen             = 'combined';                    % combine absolute and relative range scenarios             
        rangeGenType_impScen              = 'equidistant';                 % equidistant: equidistant range shifts
        scenCombType_impScen              = 'individual';                  % individual:  no combination of range and setup scenarios,
        includeNomScen_impScen            = false;                         % exclude nominal scenario
        
        % 'rndScen'   default parameters for random sampling
        numOfShiftScen_rndScen            = [25 25 25];                    % number of shifts in x y and z direction
        shiftSize_rndScen                 = [3 3 3];                       % given in [mm]
        shiftGenType_rndScen              = 'sampled';                     % sample shifts from normal distribution
        shiftCombType_rndScen             = 'combined';                    % individual:  no combination of shift scenarios;
        numOfRangeShiftScen_rndScen       = 20;                            % number of absolute and/or relative range scnearios.
        maxAbsRangeShift_rndScen          = 1;                             % maximum absolute over and undershoot in mm 
        maxRelRangeShift_rndScen          = 3.5;                           % maximum relative over and undershoot in % 
        rangeCombType_rndScen             = 'combined';                    % combine absolute and relative range scenarios             
        rangeGenType_rndScen              = 'sampled';                     % sampled: sample range shifts from normal distribution
        scenCombType_rndScen              = 'combined';                    % combine range and setup scenarios if their scenario number is consistent  
        includeNomScen_rndScen            = false;                         % exclude nominal scenario
    end
    
   %% methods
   
   % public methods go here
   methods 
       
      %Attach a listener to a TYPE 
      function attachListener(obj)
         addlistener(obj,'TYPE','PostSet',@PropLis.propChange);
      end
      
      % default constructor
      function this = matRad_multScen(ct,TYPE)
          
          if exist('TYPE','var') && ~isempty(TYPE)
              if sum(strcmp(this.AvailableScenCreationTYPE,TYPE))>0
                    this.TYPE = TYPE;
              else
                   matRad_dispToConsole(['matRad_multScen: Unknown TYPE - using the nominal scenario now'],[],'warning')
                   this.TYPE = this.DEFAULT_TYPE;
              end
          else
            this.TYPE = this.DEFAULT_TYPE;
          end
          this      = getMultScenParam(ct,this);
          this      = setMultScen(this);
          this      = calcScenProb(this);
      end % end constructor
      
      % create valid instance of an object
      function this = matRad_createValidInstance(this)
          this      = setMultScen(this);
          this      = calcScenProb(this);
      end
      
      %% setters to check for valid input
      function this = set.numOfShiftScen(this,value)
          
         if ~isempty(value)
            this.numOfShiftScen = value;
         else
            this.numOfShiftScen = this.numOfRangeShiftScen_nomScen; 
         end
         
      end
      
   end % end public methods 
    
   
   % public static methods go here, they can be called without creating an
   % instance of this class
   methods(Static)    
   
   end % end static public methods
   
   
   %private methods go here
   methods(Access = private)
       
       %%
       % set parameters according to the choosing scenario generation type
       function this = getMultScenParam(ct,this)
           
           % set all parameters according to the choosen scenario creation
           % type
           if isempty(ct)
                this.numOfCtScen = 1;   
           else
                this.numOfCtScen = ct.numOfCtScen;
           end
           
           this.numOfShiftScen       = this.(['numOfShiftScen_' this.TYPE]);
           this.shiftSize            = this.(['shiftSize_' this.TYPE]);
           this.shiftGenType         = this.(['shiftGenType_' this.TYPE]);        
           this.shiftCombType        = this.(['shiftCombType_' this.TYPE]);  
           this.numOfRangeShiftScen  = this.(['numOfRangeShiftScen_' this.TYPE]);  
           this.maxAbsRangeShift     = this.(['maxAbsRangeShift_' this.TYPE]); 
           
           this.maxRelRangeShift     = this.(['maxRelRangeShift_' this.TYPE]); 
           this.rangeCombType        = this.(['rangeCombType_' this.TYPE]); 
           this.rangeGenType         = this.(['rangeGenType_' this.TYPE]); 
           
           this.scenCombType         = this.(['scenCombType_' this.TYPE]); 
           this.includeNomScen       = this.(['includeNomScen_' this.TYPE]);  
       end
       
       
       %%
       % creates individual treatment planning scenarios
       function this = setMultScen(this)
           
            if this.includeNomScen
               nomScen = 0;
            else
               nomScen = [];
            end
            
            %% calc setup shift scenarios
            switch this.shiftGenType
                case 'equidistant'
                    % create grid vectors
                    isoShiftVec{1} = [0 linspace(-this.shiftSize(1), this.shiftSize(1), this.numOfShiftScen(1))];
                    isoShiftVec{2} = [0 linspace(-this.shiftSize(2), this.shiftSize(2), this.numOfShiftScen(2))];
                    isoShiftVec{3} = [0 linspace(-this.shiftSize(3), this.shiftSize(3), this.numOfShiftScen(3))];
                case 'sampled'
                    meanP = zeros(1,3); % mean (parameter)
                    rng('shuffle');
                    isoShiftVec{1} = [0 this.shiftSD(1) .* randn(1, this.numOfShiftScen(1)) + meanP(1)];
                    rng('shuffle');
                    isoShiftVec{2} = [0 this.shiftSD(2) .* randn(1, this.numOfShiftScen(2)) + meanP(2)];
                    rng('shuffle');
                    isoShiftVec{3} = [0 this.shiftSD(3) .* randn(1, this.numOfShiftScen(3)) + meanP(3)];     
                otherwise
                    matRad_dispToConsole('did not expect that','error');
            end

            % create scenMaskIso for isoShifts
            numIso(1) = numel(isoShiftVec{1});
            numIso(2) = numel(isoShiftVec{2});
            numIso(3) = numel(isoShiftVec{3});
            
            scenMaskIso = false(numIso(1), numIso(2), numIso(3));
            
            switch this.shiftCombType
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
                        this.shiftCombType = 'individual';
                        matRad_dispToConsole('Numnber of isoShifts in every direction has to be equal in order to perform direct combination. Performing individually instead.',[],'warning');
                        % call the function itself to get a working combination
                        [this] = setMultScen(this);
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

            % create isoShift vector based on the matrix and matching
            this.isoShift = zeros(size(matchMaskIso,1),3);
            if numel(isoShiftVec{1}) + numel(isoShiftVec{2}) + numel(isoShiftVec{3}) > 0
               for i = 1:size(matchMaskIso,1)
                   matchPos = num2cell(matchMaskIso{i,2});
                   if ~isequal([isoShiftVec{1}(matchPos{1}) isoShiftVec{2}(matchPos{2}) isoShiftVec{3}(matchPos{3})],[0 0 0]) || i == 1
                         this.isoShift(i,:) = [isoShiftVec{1}(matchPos{1}) isoShiftVec{2}(matchPos{2}) isoShiftVec{3}(matchPos{3})] * ...
                                         scenMaskIso(matchPos{:});
                   end
               end
            end
                       
            if ~this.includeNomScen
                if size(this.isoShift,1)>1
                     this.isoShift = this.isoShift(2:end,:);               % cut away the first (=nominal) scenario
                else
                     this.isoShift = [];
                end
            end
            
            this.totNumShiftScen = size(this.isoShift,1);                  % total number of shift scenarios in x,y and z direction
            
            %% calc range shift scenarios
            switch this.rangeGenType
                case 'equidistant'
                    this.relRangeShift = [nomScen linspace(-this.maxRelRangeShift, this.maxRelRangeShift, this.numOfRangeShiftScen)];
                    this.absRangeShift = [nomScen linspace(-this.maxAbsRangeShift, this.maxAbsRangeShift, this.numOfRangeShiftScen)];
                case 'sampled'
                    % relRange
                    std = this.rangeRelSD; meanP = 0;
                    rng('shuffle');
                    this.relRangeShift = [nomScen std .* randn(1, this.numOfRangeShiftScen) + meanP];
                    % absRange
                    std = this.rangeAbsSD; meanP = 0;
                    rng('shuffle');
                    this.absRangeShift = [nomScen std .* randn(1, this.numOfRangeShiftScen) + meanP];
                otherwise
                    matRad_dispToConsole('Not a valid type of generating data.','error');
            end
            
            numOfRelRangeShift = numel(this.relRangeShift);
            numOfAbsRangeShift = numel(this.absRangeShift);

            if ~isequal(numOfRelRangeShift,numOfAbsRangeShift)
                matRad_dispToConsole('Number of relative and absolute range shifts must not differ.','error');
            else
                this.totNumRangeScen = numOfRelRangeShift;
            end

            this.relRangeShift = (this.relRangeShift./100);  % consider [%]

            % check how absolute range error scenarios and relative range error scenarios should be combined
            switch this.rangeCombType    
               case 'individual'
                  rangeShift = zeros(length(this.relRangeShift)*2 ,2);
                  for i = 1:length(this.relRangeShift)
                     rangeShift((i * 2)-1,1)  = this.absRangeShift(i);
                     rangeShift((i * 2)  ,2)  = this.relRangeShift(i);
                  end

                  this.totNumRangeScen = (numOfRelRangeShift*2)-1;
                  this.absRangeShift   = [this.absRangeShift zeros(1,length(this.relRangeShift(this.relRangeShift~=0)))];
                  this.relRangeShift   = [zeros(1,length(this.relRangeShift)) this.relRangeShift(this.relRangeShift~=0)] ;

               case 'combined' 
                  rangeShift = zeros(length(this.relRangeShift),2);
                  for i = 1:1:length(this.relRangeShift) 
                     rangeShift(i,1) = this.absRangeShift(i);
                     rangeShift(i,2) = this.relRangeShift(i);
                  end

                otherwise
            end
            
            % combine setup and range scenarios according to scenCombType
            switch this.scenCombType

                case 'individual' % combine setup and range scenarios individually

                   % range errors should come first
                   if this.includeNomScen
                      this.scenForProb                                = zeros(size(this.isoShift,1)-1 + size(rangeShift,1),5);
                      this.scenForProb(1:size(rangeShift,1),4:5)      = rangeShift;
                      this.scenForProb(size(rangeShift,1)+1:end,1:3)  = this.isoShift(2:end,:);
                   else
                       this.scenForProb                               = zeros(size(this.isoShift,1) + size(rangeShift,1),5);
                       this.scenForProb(1:size(rangeShift,1),4:5)     = rangeShift;
                       this.scenForProb(size(rangeShift,1)+1:end,1:3) = this.isoShift;
                   end

                case 'permuted'

                    this.scenForProb  = zeros(size(this.isoShift,1) * size(rangeShift,1),5);
                    Cnt = 1;
                    for i = 1:size(this.isoShift,1)
                       for j = 1:size(rangeShift,1)
                          this.scenForProb(Cnt,:)   = [this.isoShift(i,:) rangeShift(j,:)];
                          Cnt = Cnt + 1;
                       end
                    end

                case 'combined'

                   if size(this.isoShift,1) == size(rangeShift,1)  && this.totNumShiftScen > 0 && this.totNumRangeScen > 0
                       this.scenForProb            = zeros(size(this.isoShift,1),5);
                       this.scenForProb(1:end,1:3) = this.isoShift;
                       this.scenForProb(1:end,4:5) = rangeShift; 
                   else
                       matRad_dispToConsole('number of setup and range scenarios MUST be the same \n',[],'warning');
                       this.scenCombType = 'individual';
                       multScen          = setMultScen(this);
                       this.scenForProb       = multScen.scenForProb;
                   end
            end
            
            % sanity check
            UniqueRowScenForProb = unique(this.scenForProb,'rows');

            if size(UniqueRowScenForProb,1) ~= size(this.scenForProb,1) && size(UniqueRowScenForProb,1)>1
                 matRad_dispToConsole('Some scenarios seem to be defined multiple times',[],'warning');
            end

            %% setup and fill combinatorics mask
            % 1st dim: ct scenarios, 
            % 2nd dim: shift scenarios, 
            % 3rd dim: range scenarios 
            this.scenMask        = false(this.numOfCtScen, this.totNumShiftScen, this.totNumRangeScen); 
            this.scenMask(:,1,1) = true; % ct scenarios

            % switch between combination modes here
            % only makes scence when numOfShiftScen>0 and numOfRangeShiftScen>0;
            if this.totNumShiftScen > 0 && this.totNumRangeScen > 0
               switch this.scenCombType
                   case 'individual'
                       % get all setup scenarios
                       [~,ixUnq] = unique(this.scenForProb(:,1:3),'rows','stable');
                       this.scenMask  = false(this.numOfCtScen, length(ixUnq), this.totNumRangeScen);
                       this.scenMask(1,:,1) = true; % iso shift scenarios
                       this.scenMask(1,1,:) = true; % range shift scenarios
                   case 'permuted'
                       this.scenMask(:,:,:) = true;
                   case 'combined'
                       % check if scenForProb is cubic (ignore ct scen) - if so fill diagonal
                       if isequal(this.totNumShiftScen, this.totNumRangeScen)
                           for i = 1:size(this.scenForProb, 1)
                               this.scenMask(1,i,i) = true;
                           end
                       else
                           this.shiftCombType = 'individual';
                           matRad_dispToConsole('number of setup and range scenarios MUST be the same \n',[],'warning');
                           this = setMultScen(this);
                       end
               end
            end

            
            % create linearalized mask where the i row points to the indexes of scenMask
            [x{1}, x{2}, x{3}] = ind2sub(size(this.scenMask),find(this.scenMask));
            this.linearMask    = cell2mat(x);
            this.totNumScen    = size(this.scenForProb,1);
            
       end % end of setMultScen
       
       
       %%
       % provides different ways of calculating the probability 
       % of occurance of individual scenarios
       function this  = calcScenProb(this)
           
            mu    = [0 0 0 0 0];
            sigma = [this.shiftSD this.rangeAbsSD this.rangeRelSD/100]; 
           
            if isequal(this.DEFAULT_probDist,'normDist')

                this.scenProb = 1;

                if isequal(this.DEFAULT_TypeProp,'probBins')

                    for i = 1:length(mu)
                        samplePosSorted = sort(unique(this.scenForProb(:,i)));
                        if numel(samplePosSorted) == 1 || sigma(i) == 0
                            continue;
                        end
                        binWidth        = (samplePosSorted(2) - samplePosSorted(1));
                        lowerBinLevel   = this.scenForProb(:,i) - 0.5*binWidth;
                        upperBinLevel   = this.scenForProb(:,i) + 0.5*binWidth;

                        this.scenProb   = this.scenProb.*0.5.*(erf((upperBinLevel-mu(i))/(sqrt(2)*sigma(i)))-erf((lowerBinLevel-mu(i))/(sqrt(2)*sigma(i))));
                    end

                elseif isequal(this.DEFAULT_TypeProp,'pointwise')
                    for i = 1:length(mu)
                        this.scenProb = this.scenProb .* (1/sqrt(2*pi*sigma(i)^2)*exp(-((this.scenForProb(:,i)-mu(i)).^2./(2*sigma(i)^2))));
                    end
                end

                % normalize probabilities since we use only a subset of
                % the 3D grid 
                this.scenProb = this.scenProb./sum(this.scenProb);

            elseif isequal(probDist,'equalProb')

               numScen  = size(samplePos,1);
               this.scenProb = repmat(1/numScen,1,numScen);

            else
                matRad_dispToConsole('Until now, only normally distributed scenarios implemented',[],'error')
            end
            
        
       end % end of calcScenProb
       
       
       
       
   end % end private methods
   
   
   
end  % end of matRad_multScen class











