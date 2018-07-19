classdef matRad_multScen < handle
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
   %                       'nomScen'    create only the nominal scenario
   %                       'wcScen'     create worst case scenarios
   %                       'impScen'    create important/grid scenarios
   %                       'rndScen'    create random scenarios
   %                       'apmScen'    create a APM scneario
   %                       'rndScenCov' creates random samples based on beam and
   %                                    covariances matrices
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
                              % 'rndScenCov' create random scenarios including ray and beam correlations
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
      
      scenMask;               % logical index array determining the activation of scenarios
      linearMask;             % linear indices of scenarios obtained from find(dij)
      
      mCovLatRnd;             % covariance matrix of random lateral erros
      mCovLatSys;             % covariance matrix of systematic lateral erros
      mCovRangeRnd;           % covariance matrix of random range erros
      mCovRangeSys;           % covariance matrix of systematic range erros
      
      vSampLatRndX;           % lateral samples from random covariance matrix
      vSampLatSysX;           % lateral samples from systematic covariance matrix
      vSampLatRndZ;           % lateral samples from random covariance matrix
      vSampLatSysZ;           % lateral samples from systematic covariance matrix
      vSampRangeRnd;          % range samples from random covariance matrix
      vSampRangeSys;          % range samples from systematic covariance matrix
      
      mCovBio;                % raywise beam correlation skelet for biological uncertainties - same as range error covariance matrix
      mCovSpot;               % component wise correlation of one individual pencil beam
      
      % values for APM
      rangeSDsys  = 3.5;      % given in %
      rangeSDrnd  = 1;        % given in [mm]
      
      shiftSDsys  = 2;        % absolut values given in [mm]
      shiftSDrnd  = 1;        % absolut values given in [mm]
      
   end
   
   % private properties which can only be changed inside matRad_multScen
   properties(SetAccess = private)
      
      DEFAULT_TYPE      = 'nomScen';
      DEFAULT_probDist  = 'normDist'; % 'normDist': normal probability distrubtion or 'equalProb' for uniform probability distribution
      DEFAULT_TypeProp  = 'probBins'; % 'probBins': cumulative prob in bin or 'pointwise' for point probability
      
      % helper objects for covariance creation
      bixelIndexRayOffset;
      bixelIndexBeamOffset;
      bixelBeamLUT;
      bixelRayLUT;
      bixelEnergyLUT;
      
   end
   
   % constant properties which are visible outside of matRad_multScen
   properties(Constant = true)
      
      AvailableScenCreationTYPE = {'nomScen','apmScen','wcScen','impScen','rndScen','rndScenCov'};
     

      corrModel   = 'block';            % either correlated,uncorrelated or block
      relBioUCT = 0.25;

      % values for random sampling and all other scenario generation types
      shiftSD     = [2.25 2.25 2.25];   % standard deviation given in [mm]
      
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
      
      % 'apmScen' default parameters for nominal scenario
      numOfShiftScen_apmScen            = [0 0 0];                       % number of shifts in x y and z direction
      shiftSize_apmScen                 = [0 0 0];                       % given in [mm]
      shiftGenType_apmScen              = 'equidistant';                 % equidistant: equidistant shifts,
      shiftCombType_apmScen             = 'individual';                  % individual:  no combination of shift scenarios;
      numOfRangeShiftScen_apmScen       = 0                              % number of absolute and/or relative range scnearios.
      maxAbsRangeShift_apmScen          = 0;                             % maximum absolute over and undershoot in mm
      maxRelRangeShift_apmScen          = 0;                             % maximum relative over and undershoot in %
      rangeCombType_apmScen             = 'combined';                    % combine absolute and relative range scenarios
      rangeGenType_apmScen              = 'equidistant';                 % equidistant: equidistant range shifts,
      scenCombType_apmScen              = 'individual';                  % individual:  no combination of range and setup scenarios,
      includeNomScen_apmScen            = true;
      
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
      
      % 'rndScen'  create random scenarios
      numOfShiftScen_rndScen            = ones(1,3) * 20;                % number of shifts in x y and z direction
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
      
      
      % 'rndScenCov'  create random scenarios
      numOfShiftScen_rndScenCov         = ones(1,3) * 200;               % number of shifts in x y and z direction
      shiftSize_rndScenCov              = [3 3 3];                       % given in [mm]
      shiftGenType_rndScenCov           = 'sampled';                     % sample shifts from normal distribution
      shiftCombType_rndScenCov          = 'combined';                    % individual:  no combination of shift scenarios;
      numOfRangeShiftScen_rndScenCov    = 200;                           % number of absolute and/or relative range scnearios.
      maxAbsRangeShift_rndScenCov       = 1;                             % maximum absolute over and undershoot in mm
      maxRelRangeShift_rndScenCov       = 3.5;                           % maximum relative over and undershoot in %
      rangeCombType_rndScenCov          = 'combined';                    % combine absolute and relative range scenarios
      rangeGenType_rndScenCov           = 'sampled';                     % sampled: sample range shifts from normal distribution
      scenCombType_rndScenCov           = 'combined';                    % combine range and setup scenarios if their scenario number is consistent
      includeNomScen_rndScenCov         = false;                         % exclude nominal scenario
   end
  
   %% methods
   
   % public methods go here
   methods
      
      %Attach a listener to a TYPE
      function attachListener(obj)
         addlistener(obj,'TYPE','PostSet',@PropLis.propChange);
      end
      
      % default constructor
      function this = matRad_multScen(ct,scenGenType)
         
         TYPE = scenGenType;
         
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
         
         this = getMultScenParam(ct,this);
         this = setMultScen(this);
         this = calcScenProb(this); 
         
      end % end constructor
      
      % create valid instance of an object
      function this = matRad_createValidInstance(this)
         this = setMultScen(this);
         this = calcScenProb(this);
      end
      
      
      function this = getCovariancesSamples(this,ct,cst,pln,stf)
        
         NumSpot    =  sum([stf(:).totalNumOfBixels]);
         
         this = setCovarianceMatrix(this,ct,cst,pln,stf);
         
         
         tmp      = full(this.mCovRangeRnd);
         mCovFull = tmp + tmp';
         mCovFull(1:NumSpot+1:end) = diag(tmp);
         [U,V,S]            = eig(mCovFull);
         this.vSampRangeRnd = bsxfun(@plus,0,(U *real(sqrtm(V))*S')' * randn(NumSpot,this.totNumScen,1))';
         
         tmp      = full(this.mCovRangeSys);
         mCovFull = tmp + tmp';
         mCovFull(1:NumSpot+1:end) = diag(tmp);
         [U,V,S]            = eig(mCovFull);
         this.vSampRangeSys = bsxfun(@plus,0,(U *real(sqrtm(V))*S')' * randn(NumSpot,this.totNumScen,1))';
         
         tmp      = full(this.mCovLatRnd);
         mCovFull = tmp + tmp';
         mCovFull(1:NumSpot+1:end) = diag(tmp);
         [U,V,S]           = eig(mCovFull);
         this.vSampLatRndX = bsxfun(@plus,0,(U *real(sqrtm(V))*S')' * randn(NumSpot,this.totNumScen,1))';
         this.vSampLatRndZ = bsxfun(@plus,0,(U *real(sqrtm(V))*S')' * randn(NumSpot,this.totNumScen,1))';
         
         tmp      = full(this.mCovLatSys);
         mCovFull = tmp + tmp';
         mCovFull(1:NumSpot+1:end) = diag(tmp);
         [U,V,S]           = eig(mCovFull);
         this.vSampLatSysX = bsxfun(@plus,0,(U *real(sqrtm(V))*S')' * randn(NumSpot,this.totNumScen,1))';
         this.vSampLatSysZ = bsxfun(@plus,0,(U *real(sqrtm(V))*S')' * randn(NumSpot,this.totNumScen,1))';
         
         this.scenProb = ones(this.totNumScen,1) * 1/this.totNumScen;

      end
      
      %
      function this = setCovarianceMatrix(this,ct,cst,pln,stf)
         
         % create APM specific error correlation matrices
         if sum(strcmp(this.TYPE,{'apmScen','rndScenCov'})) > 0
            
            options  = matRad_probOptions(ct,cst,pln);
            robID = '';
            for oo = 1:size(cst,1)
               for pp = 1:numel(cst{oo,6})
                  if ~strcmp(cst{oo,6}(pp).robustness,'none')
                     robID = 'Rob'; break;
                  end
               end
            end
            
            numComp = 10;
            
            if strcmp(pln.radiationMode,'carbon')
               numComp = 13;
            end
            
            if strcmp(robID,'') && strcmp(options.probOpt.InputUCT,'phys')
               this.mCovSpot= zeros(numComp);
            else
               for j = 1:numComp
                  for m = 1:numComp
                     this.mCovSpot(j,m) = exp(-(j-m)^2./(300));
                  end
               end
            end
            
            totalNumberOfBixels = sum([stf.totalNumOfBixels]);
            
            % obtain covariance structure
            this                                   = getBixelBeamMaps(this,pln,stf);
            [mCovlatSkelet,mCovrangeSkelet,vRange] = getCorrMatrixWrapper(this,pln,stf);
            
            % define correlations matrix for random erros
            if this.shiftSDrnd  == 0 || (strcmp(options.probOpt.InputUCT,'bio') && strcmp(robID,''))
               this.mCovLatRnd = sparse(zeros(totalNumberOfBixels));
            else
               this.mCovLatRnd = triu(mCovlatSkelet * this.shiftSDrnd^2);
            end
            
            if this.rangeSDrnd  == 0 || (strcmp(options.probOpt.InputUCT,'bio') && strcmp(robID,''))
               this.mCovRangeRnd = sparse(zeros(totalNumberOfBixels));
            else
               this.mCovRangeRnd = triu(mCovrangeSkelet * this.rangeSDrnd^2);
            end
            
            % define correlations matrix for systematic erros
            if this.shiftSDsys == 0 || (strcmp(options.probOpt.InputUCT,'bio') && strcmp(robID,''))
               this.mCovLatSys = sparse(zeros(totalNumberOfBixels));
            else
               this.mCovLatSys = triu(mCovlatSkelet * this.shiftSDsys^2);
            end

            if this.rangeSDsys  == 0 || (strcmp(options.probOpt.InputUCT,'bio') && strcmp(robID,''))
               this.mCovRangeSys = sparse(zeros(totalNumberOfBixels));
            else
               this.mCovRangeSys = triu(mCovrangeSkelet * ((this.rangeSDsys./100)^2) .* (vRange' * vRange));
            end
            
            if this.relBioUCT == 0 || (strcmp(options.probOpt.InputUCT,'phys') && strcmp(robID,''))
               this.mCovBio = sparse(zeros(totalNumberOfBixels)); 
            else
               this.mCovBio = triu(ones(totalNumberOfBixels));%mCovrangeSkelet;
            end
            
         end
      end
      
      
      
   end % end public methods
   
   
   % public static methods go here, they can be called without creating an
   % instance of this class
   methods(Static)
      
      %create correlation matrices for APM and future random sampling)
      function covMatrix = createCorrelationMatrix(stf,bixelIndexBeamOffset,bixelIndexRayOffset,errorType,corrModel)
         
         totalNumOfBixels = sum([stf.totalNumOfBixels]);
         covMatrix        = spalloc(totalNumOfBixels,totalNumOfBixels,1);
         
         switch corrModel
            case {'correlated'}
               
               warning('creating correlated covariance matrixes may lead to memory overflow');
               %covMatrix = spalloc(pln.uct.totalNumOfBixels,pln.uct.totalNumOfBixels,pln.uct.totalNumOfBixels*pln.uct.totalNumOfBixels)
               covMatrix  = ones(totalNumOfBixels);
               
            case {'block'}
               
               if strcmp(errorType,'setup')
                  % bixels belonging to the same beam direction are perfectly correlated
                  for idx = 2:length(bixelIndexBeamOffset)
                     covMatrix(bixelIndexBeamOffset(idx-1):bixelIndexBeamOffset(idx),...
                        bixelIndexBeamOffset(idx-1):bixelIndexBeamOffset(idx)) = 1; %#ok<*SPRIX>
                  end
                  
               elseif strcmp(errorType,'range')
                  % bixels belonging to the same beam and impinging at the
                  % same lateral position having different energies are
                  % considered perfectly correlated as they penetrate the
                  % same tissue on their ray
                  for idx = 2:1:length(bixelIndexRayOffset)
                     
                     if idx == 2
                        covMatrix(bixelIndexRayOffset(idx-1):bixelIndexRayOffset(idx),...
                           bixelIndexRayOffset(idx-1):bixelIndexRayOffset(idx)) = 1;
                     else
                        covMatrix(bixelIndexRayOffset(idx-1)+1:bixelIndexRayOffset(idx),...
                           bixelIndexRayOffset(idx-1)+1:bixelIndexRayOffset(idx)) = 1;
                     end
                     
                  end
                  
               end
               
            case {'uncorrelated'}
               covMatrix = speye(totalNumOfBixels,totalNumOfBixels);
         end
         
      end %eof createCorrelationMatrix
      
      
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
               std = this.rangeSDsys; meanP = 0;
               rng('shuffle');
               this.relRangeShift = [nomScen std .* randn(1, this.numOfRangeShiftScen) + meanP];
               % absRange
               std = this.rangeSDrnd; meanP = 0;
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
         sigma = [this.shiftSD this.rangeSDrnd this.rangeSDsys/100];
         
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
      
      
      %%
      % get helper objects for the creation of uncertainty covariance matrices
      function this = getBixelBeamMaps(this,pln,stf)
         
         % calculate total number of bixels from steering file
         totalNumOfBixels     = sum([stf(:).totalNumOfBixels]);
         
         %Get for each bixel the beam index and ray index
         this.bixelIndexBeamOffset = [1 cumsum([stf(:).totalNumOfBixels])];
         this.bixelIndexRayOffset  = [1 cumsum([stf(:).numOfBixelsPerRay])];
         
         % create spot beam and spot ray look up tables
         this.bixelBeamLUT   = zeros(totalNumOfBixels,1);
         this.bixelRayLUT    = zeros(totalNumOfBixels,1);
         this.bixelEnergyLUT = zeros(totalNumOfBixels,1);
         
         cntBeam = 1;  cntRay = 1; cntBixel = 1;
         % loop over all beams
         for i = 1:pln.propStf.numOfBeams
            this.bixelBeamLUT(cntBeam:cntBeam+stf(i).totalNumOfBixels-1) = i;
            cntBeam = cntBeam + stf(i).totalNumOfBixels;
            % loop over all rays per beam
            for j = 1:stf(i).numOfRays
               numSpots = length(stf(i).ray(j).energy);
               this.bixelRayLUT(cntBixel:cntBixel + numSpots -1 ) = cntRay;
               this.bixelEnergyLUT(cntBixel:cntBixel + numSpots -1) = stf(i).ray(j).energyIx;
               cntRay = cntRay + 1;
               cntBixel = cntBixel + numSpots;
            end
            cntRay = 1;
         end
         
         
      end %eof getBixelBeamMaps
      
      %%
      % calculate covariance matrices for setup and range error
      function [mCovlatSkelet,mCovrangeSkelet,vRange] = getCorrMatrixWrapper(this,pln,stf)
         
         mCovlatSkelet   =  this.createCorrelationMatrix(stf,this.bixelIndexBeamOffset,this.bixelIndexRayOffset,'setup',this.corrModel) ;
         mCovrangeSkelet =  this.createCorrelationMatrix(stf,this.bixelIndexBeamOffset,this.bixelIndexRayOffset,'range',this.corrModel) ;
         
         % obtain ranges
         % prepare structures necessary for particles
         fileName = [pln.radiationMode '_' pln.machine];
         try
            load([fileparts(mfilename('fullpath')) filesep fileName]);
         catch
            matRad_dispToConsole(['Could not find the following machine file: ' fileName ],param,'error');
         end
         
         vRange = [];               % range of bixel in mm
         for i = 1:size(stf,2)
            for j = 1:stf(i).numOfRays
               vRange = [vRange round([machine.data(stf(i).ray(j).energyIx).range])];  %#ok<AGROW>
            end
         end
         
         %figure,
         %subplot(121),spy(this.mCOVlat)
         %subplot(121),spy(this.mCOVrange)
         
      end %eof
      
      
   end % end private methods
   
   
   
end  % end of matRad_multScen class











