classdef matRad_multScen
    %  matRad_multScen
    %  This class creates all required biological model parameters according to
    % a given radiation modatlity and a given bio model identifier string.
    %
    % constructor call
    %   pln.multScen = matRad_multScen(ct,TYPE,uctModel,combineScenarios);
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
        
        TYPE;                   % denotes the sampling type which can be one of the following
        % 'nomScen'   create only the nominal scenario
        % 'wcScen'    create worst case scenarios
        % 'impScen'   create important/grid scenarios
        % 'rndScen'   create random scenarios
        
        
        probDist;               % denotes probability distrubtion which can be one of the following
        % 'normDist': normal probability distrubtion
        % 'equalProb' for uniform probability distribution
        
        typeProp;
        
        % Uncertainty Model
        rangeRelSD  = 3.5;                % given in %
        rangeAbsSD  = 1;                  % given in [mm]
        shiftSD     = [2.25 2.25 2.25];   % given in [mm]
        wcFactor    = 1.5;                % Multiplier to compute the worst case / maximum shifts
        
        
        shiftGenType;           % 'equidistant': equidistant shifts, 'sampled': sample shifts from normal distribution
        shiftCombType;          % 'individual':  no combination of shift scenarios;       number of shift scenarios is sum(multScen.numOfShiftScen)
        % 'combined':    combine shift scenarios;                 number of shift scenarios is numOfShiftScen
        % 'permuted':    create every possible shift combination; number of shift scenarios is 8,27,64 .
        
        rangeCombType;          % 'individual':  no combination of absolute and relative range scenarios
        % 'combined':    combine absolute and relative range scenarios
        rangeGenType;           % 'equidistant': equidistant range shifts, 'sampled': sample range shifts from normal distribution
        
        scenCombType;           % 'individual':  no combination of range and setup scenarios,
        % 'combined':    combine range and setup scenarios if their scenario number is consistent
        % 'permuted':    create every possible combination of range and setup scenarios
        
        includeNomScen;         % boolean to determine if the nominal scenario should be included
        
        numOfShiftScen;         % number of shifts in x y and z direction
        numOfRangeShiftScen;    % number of absolute and/or relative range scnearios.
    end
    
    %These properties are accessible but not editable
    properties (SetAccess = private)
        
        % ct scenarios
        numOfCtScen = 1;            % number of imported ct scenarios
        
        % shift scenarios
        
        shiftSize;              % 3x1 vector to define maximal shift in [mm]  % (e.g. abdominal cases 5mm otherwise 3mm)
        
        % b) define range error scenarios
        
        % if absolute and relative range scenarios are defined then rangeCombType defines the resulting number of range scenarios
        maxAbsRangeShift;       % maximum absolute over and undershoot in mm
        maxRelRangeShift;       % maximum relative over and undershoot in %
        
        
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
    end
    
    % constant properties which are visible outside of matRad_multScen
    properties(Constant = true)
        
        AvailableScenCreationTYPE = {'nomScen','wcScen','impScen','rndScen'};
        AvailableProbDist = {'normDist','equalProb'};
        AvailableTypeProp = {'probBins','pointwise'};
        
    end
    
    % constant deafult properties which are only visible within matRad_multScen
    properties (Constant = true, Access = private)
        DEFAULT_TYPE      = 'nomScen';
        DEFAULT_probDist  = 'normDist'; % 'normDist': normal probability distrubtion or 'equalProb' for uniform probability distribution
        DEFAULT_TypeProp  = 'probBins'; % 'probBins': cumulative prob in bin or 'pointwise' for point probability
        
        % 'nomScen' default parameters for nominal scenario
        numOfShiftScen_nomScen            = [0 0 0];                       % number of shifts in x y and z direction
        %shiftSize_nomScen                 = [0 0 0];                       % given in [mm]
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
        %shiftSize_wcScen                  = [4 4 4];                       % given in [mm]
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
        numOfShiftScen_impScen            = [10 10 10];                       % number of shifts in x y and z direction
        %shiftSize_impScen                 = [0 0 0];                       % given in [mm]
        shiftGenType_impScen              = 'equidistant';                 % equidistant: equidistant shifts
        shiftCombType_impScen             = 'individual';                  % individual:  no combination of shift scenarios
        numOfRangeShiftScen_impScen       = 10;                            % number of absolute and/or relative range scnearios
        maxAbsRangeShift_impScen          = 1;                             % maximum absolute over and undershoot in mm
        maxRelRangeShift_impScen          = 3.5;                           % maximum relative over and undershoot in %
        rangeCombType_impScen             = 'combined';                    % combine absolute and relative range scenarios
        rangeGenType_impScen              = 'equidistant';                 % equidistant: equidistant range shifts
        scenCombType_impScen              = 'individual';                  % individual:  no combination of range and setup scenarios,
        includeNomScen_impScen            = false;                         % exclude nominal scenario
        
        % 'rndScen'   default parameters for random sampling
        numOfShiftScen_rndScen            = [25 25 25];                    % number of shifts in x y and z direction
        %shiftSize_rndScen                 = [3 3 3];                       % given in [mm]
        shiftGenType_rndScen              = 'sampled';                     % sample shifts from normal distribution
        shiftCombType_rndScen             = 'combined';                    % individual:  no combination of shift scenarios;
        numOfRangeShiftScen_rndScen       = 25;                            % number of absolute and/or relative range scnearios.
        maxAbsRangeShift_rndScen          = NaN;                             % maximum absolute over and undershoot in mm
        maxRelRangeShift_rndScen          = NaN;                           % maximum relative over and undershoot in %
        rangeCombType_rndScen             = 'combined';                    % combine absolute and relative range scenarios
        rangeGenType_rndScen              = 'sampled';                     % sampled: sample range shifts from normal distribution
        scenCombType_rndScen              = 'combined';                    % combine range and setup scenarios if their scenario number is consistent
        includeNomScen_rndScen            = false;                         % exclude nominal scenario
    end
    
    properties (Access = private)
        lockUpdate = false;
        lockInit = false;
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
            
            matRad_cfg = MatRad_Config.instance();
            
            if isempty(ct)
                this.numOfCtScen = 1;
            else
                this.numOfCtScen = ct.numOfCtScen;
            end
            
            this.lockInit = true;
            if exist('TYPE','var') && ~isempty(TYPE)
                if sum(strcmp(this.AvailableScenCreationTYPE,TYPE))>0
                    this.TYPE = TYPE;
                else
                    matRad_cfg.dispWarning('matRad_multScen: Unknown TYPE - using the nominal scenario now');
                    this.TYPE = this.DEFAULT_TYPE;
                end
            else
                this.TYPE = this.DEFAULT_TYPE;
            end
            
            if exist('probDist','var') && ~isempty(probDist)
                if sum(strcmp(this.AvailableProbDist,probDist))>0
                    this.probDist = probDist;
                else
                    matRad_cfg.dispWarning('matRad_multScen: Unknown probability distribution - using normal distribution now');
                    this.probDist = this.DEFAULT_probDist;
                end
            else
                this.probDist = this.DEFAULT_probDist;
            end
            
            if exist('typeProp','var') && ~isempty(typeProp)
                if sum(strcmp(this.AvailableTypeProp,typeProp))>0
                    this.typeProp = typeProp;
                else
                    matRad_cfg.dispWarning('matRad_multScen: Unknown type prop - using the probBins now');
                    this.typeProp = this.DEFAULT_TypeProp;
                end
            else
                this.typeProp = this.DEFAULT_TypeProp;
            end
            
            this.lockInit = false;
            
            this      = this.initialize();
            
        end % end constructor
        
        function newInstance = extractSingleNomScen(this,ctScen,i)
            newInstance = this;
            newInstance.TYPE = 'nomScen';
            newInstance.lockInit = true;
            newInstance.lockUpdate = true;
            
            newInstance.scenForProb         = this.scenForProb(i,:);
            newInstance.relRangeShift       = this.scenForProb(i,5);
            newInstance.absRangeShift       = this.scenForProb(i,4);
            newInstance.isoShift            = this.scenForProb(i,1:3);
            newInstance.totNumShiftScen     = 1;
            newInstance.totNumRangeScen     = 1;
            newInstance.numOfCtScen         = ctScen;
            newInstance.scenMask            = 1;
            newInstance.linearMask          = 1;
            newInstance.scenProb            = 1;
            
            newInstance.lockInit = false;
            newInstance.lockUpdate = false;
        end
        
        % create valid instance of an object
        function this = matRad_createValidInstance(this)
            this      = setMultScen(this);
            this      = calcScenProb(this);
        end
        
        %% Setter functions
        function this = set.TYPE(this,value)
            if ischar(value) && any(strcmp(value,this.AvailableScenCreationTYPE))
                this.TYPE = value;
                this = this.initialize();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid scenario TYPE!');
            end
        end
        
        function this = set.probDist(this,value)
            if ischar(value) && any(strcmp(value,this.AvailableProbDist))
                this.probDist = value;
                this = this.initialize();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid probability distribution!');
            end
        end
         
        function this = set.typeProp(this,value)
            if ischar(value) && any(strcmp(value,this.AvailableTypeProp))
                this.typeProp = value;
                this = this.initialize();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid scenario type prop!');
            end
        end
        
        function this = set.rangeRelSD(this,value)
            if isnumeric(value) && isscalar(value) && value >= 0
                this.rangeRelSD = value;
                this = this.updateScenariosFromUncertainties();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for relative range uncertainty!');
            end
        end
        
        function this = set.rangeAbsSD(this,value)
            if isnumeric(value) && isscalar(value) && value >= 0
                this.rangeAbsSD = value;
                this = this.updateScenariosFromUncertainties();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for absolute range uncertainty!');
            end
        end
        
        function this = set.shiftSD(this,value)
            if isnumeric(value) && isvector(value) && numel(value) == 3 && all(value >= 0)
                this.shiftSD = value;
                this = this.updateScenariosFromUncertainties();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for setup uncertainty!');
            end
        end
        
        function this = set.wcFactor(this,value)
            if isnumeric(value) && isscalar(value) && value >= 0
                this.wcFactor = value;
                this = this.updateScenariosFromUncertainties();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for worst case factor!');
            end
        end
        
        function this = set.numOfShiftScen(this,value)
            if isnumeric(value) && isvector(value) && numel(value) == 3 && all(value >= 0) && all(mod(value,1) == 0)
                this.numOfShiftScen = value;
                this = this.updateScenariosFromUncertainties();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for shift scenario number!');
            end
        end
        
        function this = set.numOfRangeShiftScen(this,value)
            if isnumeric(value) && isscalar(value) && value >= 0 && mod(value,1) == 0
                this.numOfRangeShiftScen = value;
                this = this.updateScenariosFromUncertainties();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for number of range over/undershoots!');
            end
        end
        
        function this = set.shiftGenType(this,value)
            if ischar(value) && any(strcmp(value,{'equidistant','sampled','sampled_truncated'}))
                this.shiftGenType = value;
                this = this.updateScenariosFromUncertainties();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for shift generation type!');
            end
        end
        
        function this = set.rangeGenType(this,value)
            if ischar(value) && any(strcmp(value,{'equidistant','sampled','sampled_truncated'}))
                this.rangeGenType = value;
                this = this.updateScenariosFromUncertainties();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for range shift generation type!');
            end
        end
        
        function this = set.shiftCombType(this,value)
            if ischar(value) && any(strcmp(value,{'individual','combined','permuted','permuted_truncated'}))
                this.shiftCombType = value;
                this = this.updateScenariosFromUncertainties();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for shift combination!');
            end
        end
        
        function this = set.scenCombType(this,value)
            if ischar(value) && any(strcmp(value,{'individual','combined','permuted'}))
                this.scenCombType = value;
                this = this.updateScenariosFromUncertainties();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for scenario combination!');
            end
        end
        
        function this = set.rangeCombType(this,value)
            if ischar(value) && any(strcmp(value,{'individual','combined'}))
                this.rangeCombType = value;
                this = this.updateScenariosFromUncertainties();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for absolute & relative range combination!');
            end
        end
        
        function this = set.includeNomScen(this,value)
            if isnumeric(value) || islogical(value)
                this.includeNomScen = logical(value);
                this = this.updateScenariosFromUncertainties();
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Invalid value for inclusion of nominal scenario!');
            end
        end
        
        function listAllScenarios(this)
            matRad_cfg = MatRad_Config.instance();
            matRad_cfg.dispInfo('Listing all scenarios...\n');
            matRad_cfg.dispInfo('\t#\txShift\tyShift\tzShift\tabsRng\trelRng\tprob.\n');
            for s = 1:size(this.scenForProb,1)
                str = num2str(this.scenForProb(s,:),'\t%.3f');
                matRad_cfg.dispInfo('\t%d\t%s\t%.3f\n',s,str,this.scenProb(s));
            end
        end
    end % end public methods
    
    %private methods go here
    methods(Access = private)
        %%
        % set parameters according to the choosing scenario generation type
        function this = initialize(this)
            if this.lockInit
                return;
            end
            
            
            this.lockUpdate = true;
            
            this.numOfShiftScen       = this.(['numOfShiftScen_' this.TYPE]);
            this.shiftGenType         = this.(['shiftGenType_' this.TYPE]);
            this.shiftCombType        = this.(['shiftCombType_' this.TYPE]);
            this.numOfRangeShiftScen  = this.(['numOfRangeShiftScen_' this.TYPE]);
            
            this.rangeCombType        = this.(['rangeCombType_' this.TYPE]);
            this.rangeGenType         = this.(['rangeGenType_' this.TYPE]);
            
            this.scenCombType         = this.(['scenCombType_' this.TYPE]);
            this.includeNomScen       = this.(['includeNomScen_' this.TYPE]);
            
            this.lockUpdate = false;
            
            this = this.updateScenariosFromUncertainties();
        end
        
        %%
        function this = updateScenariosFromUncertainties(this)
            if this.lockUpdate
                return;
            end
            switch this.TYPE
                case 'wcScen'
                    this.shiftSize = this.shiftSD * this.wcFactor;
                    this.maxAbsRangeShift = this.rangeAbsSD * this.wcFactor;
                    this.maxRelRangeShift = this.rangeRelSD * this.wcFactor;
                case 'nomScen'
                    this.shiftSize = [0 0 0];
                    this.maxAbsRangeShift = 0;
                    this.maxRelRangeShift = 0;
                case 'impScen'
                    this.shiftSize = this.shiftSD * this.wcFactor;
                    this.maxAbsRangeShift = this.rangeAbsSD * this.wcFactor;
                    this.maxRelRangeShift = this.rangeRelSD * this.wcFactor;
                case 'rndScen'
                    this.shiftSize = [NaN NaN NaN];
                    this.maxAbsRangeShift = this.rangeAbsSD * 3;
                    this.maxRelRangeShift = this.rangeRelSD * 3;
                    
            end
            
            this = this.setMultScen();
            this = this.calcScenProb();
        end
        
        %%
        % creates individual treatment planning scenarios
        function this = setMultScen(this)
            matRad_cfg = MatRad_Config.instance();
            if this.includeNomScen
                nomScen = 0;
            else
                nomScen = [];
            end
            
            %% calc setup shift scenarios
            switch this.shiftGenType
                case 'equidistant'
                    % create grid vectors
                    if(this.numOfShiftScen(1)>0)
                        xVec = linspace(0, this.shiftSize(1), (this.numOfShiftScen(1)+1)/2);
                        xVec(1)=[];
                        isoShiftVec{1} = [0 -xVec xVec];
                    else
                        isoShiftVec{1}=0;
                    end
                    if(this.numOfShiftScen(2)>0)
                        yVec = linspace(0, this.shiftSize(2), (this.numOfShiftScen(2)+1)/2);
                        yVec(1)=[];
                        isoShiftVec{2} = [0 -yVec yVec];
                    else
                        isoShiftVec{2}=0;
                    end
                    if(this.numOfShiftScen(3)>0)
                        zVec = linspace(0, this.shiftSize(3), (this.numOfShiftScen(3)+1)/2);
                        zVec(1)=[];
                        isoShiftVec{3} = [0 -zVec zVec];
                    else
                        isoShiftVec{3}=0;
                    end
                case 'sampled'
                    meanP = zeros(1,3); % mean (parameter)
                    if matRad_getEnvironment == 'MATLAB' rng('shuffle'), end;
                    isoShiftVec{1} = [0 this.shiftSD(1) .* randn(1, this.numOfShiftScen(1)) + meanP(1)];
                    if matRad_getEnvironment == 'MATLAB' rng('shuffle'), end;
                    isoShiftVec{2} = [0 this.shiftSD(2) .* randn(1, this.numOfShiftScen(2)) + meanP(2)];
                    if matRad_getEnvironment == 'MATLAB' rng('shuffle'), end;
                    isoShiftVec{3} = [0 this.shiftSD(3) .* randn(1, this.numOfShiftScen(3)) + meanP(3)];
                case 'sampled_truncated'
                    meanP = zeros(1,3); % mean (parameter)
                    if matRad_getEnvironment == 'MATLAB' rng('shuffle'), end;
                    pd{1} = makedist('Normal','mu',meanP(1),'sigma',this.shiftSD(1));
                    t{1} = truncate(pd{1},-this.shiftSD(1) * this.wcFactor,this.shiftSD(1) * this.wcFactor);
                    isoShiftVec{1} = [0 transpose(random(t{1},this.numOfShiftScen(1),1))];
                    if matRad_getEnvironment == 'MATLAB' rng('shuffle'), end;
                    pd{2} = makedist('Normal','mu',meanP(2),'sigma',this.shiftSD(2));
                    t{2} = truncate(pd{2},-this.shiftSD(2) * this.wcFactor,this.shiftSD(2) * this.wcFactor);
                    isoShiftVec{2} = [0 transpose(random(t{2},this.numOfShiftScen(2),1))];
                    if matRad_getEnvironment == 'MATLAB' rng('shuffle'), end;
                    pd{3} = makedist('Normal','mu',meanP(3),'sigma',this.shiftSD(3));
                    t{3} = truncate(pd{3},-this.shiftSD(3) * this.wcFactor,this.shiftSD(3) * this.wcFactor);
                    isoShiftVec{3} = [0 transpose(random(t{3},this.numOfShiftScen(3),1))];
                    
                otherwise
                    matRad_cfg.dispError('did not expect that!');
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
                case 'permuted_truncated'
                    for ix=1:numel(isoShiftVec{1})
                        for iy=1:numel(isoShiftVec{2})
                            for iz=1:numel(isoShiftVec{3})
                                if isoShiftVec{1}(ix)^2/this.shiftSD(1)^2+isoShiftVec{2}(iy)^2/this.shiftSD(2)^2+isoShiftVec{3}(iz)^2/this.shiftSD(3)^2<=this.wcFactor^2;
                                    scenMaskIso(ix,iy,iz) = true;
                                end
                            end
                        end
                    end
                case 'combined'
                    % determine that matrix is cubic
                    if isequal(numIso(1), numIso(2), numIso(3))
                        for i = 1:numIso(1)
                            scenMaskIso(i,i,i) = true;
                        end
                    else
                        this.shiftCombType = 'individual';
                        matRad_cfg.dispWarning('Numnber of isoShifts in every direction has to be equal in order to perform direct combination. Performing individually instead.');
                        % call the function itself to get a working combination
                        [this] = setMultScen(this);
                    end
                otherwise
                    matRad_cfg.dispError('Uncaught exception. Probably TYPO.','error');
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
                    if matRad_getEnvironment == 'MATLAB' rng('shuffle'), end;
                    this.relRangeShift = [nomScen std .* randn(1, this.numOfRangeShiftScen) + meanP];
                    % absRange
                    std = this.rangeAbsSD; meanP = 0;
                    if matRad_getEnvironment == 'MATLAB' rng('shuffle'), end;
                    this.absRangeShift = [nomScen std .* randn(1, this.numOfRangeShiftScen) + meanP];
                case 'sampled_truncated'
                    std = this.rangeRelSD; meanP = 0;
                    % relRange
                    if matRad_getEnvironment == 'MATLAB' rng('shuffle'), end;
                    pd{1} = makedist('Normal','mu',meanP(1),'sigma',std);
                    t{1} = truncate(pd{1},-std * this.wcFactor,std * this.wcFactor);
                    this.relRangeShift = [nomScen random(t{1},this.numOfShiftScen(1),1)];
                    % absRange
                    if matRad_getEnvironment == 'MATLAB' rng('shuffle'), end;
                    pd{2} = makedist('Normal','mu',meanP(1),'sigma',std);
                    t{2} = truncate(pd{2},-std * this.wcFactor,std * this.wcFactor);
                    this.absRangeShift = [nomScen random(t{1},this.numOfShiftScen(1),1)];
                otherwise
                    matRad_cfg.dispError('Not a valid type of generating data.');
            end
            
            numOfRelRangeShift = numel(this.relRangeShift);
            numOfAbsRangeShift = numel(this.absRangeShift);
            
            if ~isequal(numOfRelRangeShift,numOfAbsRangeShift)
                matRad_cfg.dispError('Number of relative and absolute range shifts must not differ.');
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
                    
                    nIso = size(this.isoShift,1);
                    nRange = size(rangeShift,1);
                    
                    
                    % range errors should come last
                    if this.includeNomScen
                        nIso = nIso - 1;
                        nRange = nRange -1;
                        
                        this.scenForProb = zeros(1,5);  %Nominal Scenario
                        
                        rangeScen = zeros(nRange,5);
                        rangeScen(:,4:5) = rangeShift(2:end,:);
                        
                        isoScen = zeros(nIso,5);
                        isoScen(:,1:3) = this.isoShift(2:end,:);
                        this.scenForProb = [this.scenForProb; isoScen; rangeScen];
                    else
                        rangeScen = zeros(nRange,5);
                        isoScen = zeros(nIso,5);
                        rangeScen(:,4:5) = rangeShift;
                        isoScen(:,1:3) = this.isoShift;
                        
                        this.scenForProb = [isoScen; rangeScen];
                    end
                    
                case 'permuted'
                    
                    this.scenForProb  = zeros(size(this.isoShift,1) * size(rangeShift,1),5);
                    Cnt = 1;
                    for j = 1:size(rangeShift,1)
                        for i = 1:size(this.isoShift,1)
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
                        matRad_cfg.dispWarning('number of setup and range scenarios MUST be the same \n');
                        this.scenCombType = 'individual';
                        multScen          = setMultScen(this);
                        this.scenForProb       = multScen.scenForProb;
                    end
            end
            
            % sanity check
            UniqueRowScenForProb = unique(this.scenForProb,'rows');
            
            if size(UniqueRowScenForProb,1) ~= size(this.scenForProb,1) && size(UniqueRowScenForProb,1)>1
                matRad_cfg.dispWarning('Some scenarios seem to be defined multiple times');
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
                        %                        [~,ixUnq] = unique(this.scenForProb(:,1:3),'rows','stable');
                        uq   = this.scenForProb(1,1:3);
                        ixUnq = [1];
                        for col  = 2: size(this.scenForProb(:,1:3),1)
                            flag=1;
                            for colq = 1: size(uq,1)
                                if (sum(this.scenForProb(col,1:3) ==uq(colq,:)) == 3 )
                                    flag=0;
                                end
                            end
                            if flag
                                ixUnq = [ixUnq; col];
                            end
                        end
                        this.scenMask  = false(this.numOfCtScen, length(ixUnq), this.totNumRangeScen);
                        this.scenMask(:,1,1) = true; % ct scenarios
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
                            matRad_cfg.dispWarning('number of setup and range scenarios MUST be the same \n');
                            this = setMultScen(this);
                        end
                end
            end
            
            
            % create linearalized mask where the i row points to the indexes of scenMask
            [x{1}, x{2}, x{3}] = ind2sub(size(this.scenMask),find(this.scenMask));
            this.linearMask    = cell2mat(x);
            this.totNumScen    = sum(this.scenMask(:));
            
        end % end of setMultScen
        
        
        %%
        % provides different ways of calculating the probability
        % of occurance of individual scenarios
        function this  = calcScenProb(this)
            
            mu    = [0 0 0 0 0];
            sigma = [this.shiftSD this.rangeAbsSD this.rangeRelSD/100];
            
            if isequal(this.probDist,'normDist')
                
                this.scenProb = 1;
                
                if isequal(this.typeProp,'probBins')
                    
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
                    
                elseif isequal(this.typeProp,'pointwise')
                    for i = 1:length(mu)
                        this.scenProb = this.scenProb .* (1/sqrt(2*pi*sigma(i)^2)*exp(-((this.scenForProb(:,i)-mu(i)).^2./(2*sigma(i)^2))));
                    end
                end
                
                % normalize probabilities since we use only a subset of
                % the 3D grid
                this.scenProb = this.scenProb./sum(this.scenProb);
                
            elseif isequal(this.probDist,'equalProb')
                this.scenProb = repmat(1/this.totNumScen,1,this.totNumScen);
                
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('Until now, only normally distributed scenarios implemented');
            end
            
            
        end % end of calcScenProb
        
        
        
        
    end % end private methods
    
    
    
end  % end of matRad_multScen class










