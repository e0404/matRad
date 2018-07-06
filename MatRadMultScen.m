classdef MatRadMultScen < handle 
  properties               
           
      % scenarios
      numOfCtScen;            % number of imported ct scenarios                        
      % shift scenarios
      shiftSize;              % 3x1 vector to define maximal shift in [mm]  % (e.g. abdominal cases 5mm otherwise 3mm)
      maxAbsRangeShift;       % maximum absolute over and undershoot in mm   
      maxRelRangeShift;       % maximum relative over and undershoot in % 
      
      
      % combination
      rangeCombType;          % 'permuted':    create every possible range combination
                              % 'combined':    combine absolute and relative range scenario
      shiftCombType;          % 'permuted':    create every possible shift combination; 
      scenCombType;           % 'permuted':    create every possible combination of range and setup scenarios
                                               
      samplingMethod;           % 'random' / 'grid'
      includeNomScen;         % boolean to determine if the nominal scenario should be included

      totNumScen;             % total number of scenarios (might be changed later)
      numOfFrac;              % number of fractions
      
      shiftSD_sys     = [0 0 0];   % given in [mm]
      shiftSD_rand    = [2 2 2];   % given in [mm]
      
      rangeRelSD_sys  = 1.8;       % given in %
      rangeRelSD_rand = 0;         % given in %
      rangeAbsSD_sys  = 0.8;       % given in [mm]
      rangeAbsSD_rand = 0;         % given in [mm]

  end
  
  % private properties which can only be changed inside matRad_multScen
  properties (SetAccess = private)
      fraction_index % determines which dose j belongs to which treatment scenario i
      scenForProb;            % matrix for probability calculation - each row denotes one scenario  
      % these parameters will be filled according to the choosen scenario type
      isoShift;
      relRangeShift;
      absRangeShift;
      scenMask;
      
      DEFAULT_TYPE      = 'nomScen'; 
      DEFAULT_probDist  = 'normDist'; % 'normDist': normal probability distrubtion or 'equalProb' for uniform probability distribution  
      DEFAULT_TypeProp  = 'probBins'; % 'probBins': cumulative prob in bin or 'pointwise' for point probability  
  end
  
  % constant properties which are visible outside of matRad_multScen
  properties (Constant = true)
      AvailableScenCreationTYPE = {'nomScen','wcScen','impScen','rndScen'};
  end
  
  properties (Dependent = true)
      SD_sys_all;
      SD_rand_all;
      scenProb;               % probability of each scenario stored in a vector
      numberOfTreatmentScen;
  end
  
  methods
    function this = MatRadMultScen()
        
    end
    
    function randomSampling(this, totNumScen, numOfFrac, shiftSD_sys, shiftSD_rand, rangeRelSD_sys, rangeRelSD_rand, rangeAbsSD_sys, rangeAbsSD_rand)
        this.totNumScen = totNumScen;
        this.numOfFrac = numOfFrac;
        this.shiftSD_sys = shiftSD_sys;
        this.shiftSD_rand = shiftSD_rand;
        this.rangeRelSD_sys = rangeRelSD_sys;
        this.rangeRelSD_rand = rangeRelSD_rand;
        this.rangeAbsSD_sys = rangeAbsSD_sys;
        this.rangeAbsSD_rand = rangeAbsSD_rand;
        
        this.createRandomSampledScenarios();
    end
    
    function w = get.scenProb(this)
        w = this.computeScenProb();
    end
    
    function numberOfTreatmentScen = get.numberOfTreatmentScen(this)
        numberOfTreatmentScen = ceil(this.totNumScen / this.numOfFrac);
    end
    function SD_sys_all = get.SD_sys_all(this)
        SD_sys_all = [this.shiftSD_sys this.rangeAbsSD_sys this.rangeRelSD_sys];
    end
    
    function SD_rand_all = get.SD_rand_all(this)
        SD_rand_all = [this.shiftSD_rand this.rangeAbsSD_rand this.rangeRelSD_rand];
    end
    
  end % eof public methods
  
  methods (Access = private)

    function createRandomSampledScenarios(this)
        this.totNumScen = this.numberOfTreatmentScen * this.numOfFrac;
        
        this.scenForProb = NaN * ones(this.totNumScen, numel(this.SD_sys_all));
        this.fraction_index = NaN * ones(this.totNumScen, 1);

        runIx = 1;
        for i = 1:this.numberOfTreatmentScen
          mu = this.sampleFromNormalDist(0, this.SD_sys_all, size(this.SD_sys_all));
          for j = 1:this.numOfFrac
            this.scenForProb(runIx, :) = this.sampleFromNormalDist(mu, this.SD_rand_all, size(this.SD_rand_all));
            this.fraction_index(runIx) = i;
            runIx = runIx + 1;
          end
        end
        this.isoShift = this.scenForProb(:,1:3);
        this.absRangeShift = this.scenForProb(:,4);
        this.relRangeShift = this.scenForProb(:,5);
        this.scenMask = ones(this.numOfCtScen, size(this.isoShift, 1), numel(this.absRangeShift));
    end
    
    function scenProb = computeScenProb(this)
        if strcmp(this.samplingMethod, 'random')
            scenProb = 1 / this.totNumScen * ones(this.totNumScen, 1);
        else
            error('Do that later');
        end
    end
    
  end % eof private methods
  
  methods (Static)
    function v = sampleFromNormalDist(mu,sigma, dim)
      v = sigma .* randn(dim(1),dim(2)) + mu;
    end
  end
  
end
