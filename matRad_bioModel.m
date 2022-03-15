classdef matRad_bioModel
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  matRad_bioModel
    %  This class creates all required biological model parameters according to
    % a given radiation modatlity and a given bio model identifier string.
    %
    % constructor call
    %   pln.bioParam = matRad_bioModel(sRadiationMode,sQuantityOpt, sModel)
    %
    %   e.g. pln.bioParam = matRad_bioModel('protons','constRBE_RBExD')
    %
    % input
    %   sRadiationMode:     radiation modality 'photons' 'protons' 'carbon'
    %   sQuntityOpt:        string to denote the quantity used for
    %                       optimization 'physicalDose', 'RBExD', 'effect'
    %   sModel:             string to denote which biological model is used
    %                       'none': for photons, protons, carbon                                    'constRBE': constant RBE for photons and protons 
    %                       'MCN': McNamara-variable RBE model for protons                          'WED': Wedenberg-variable RBE model for protons 
    %                       'LEM': Local Effect Model for carbon ions
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
    
    % public properties which can be changed outside this class
    properties
        
    end
    
    % public properties which can only be changed inside this class
    properties(SetAccess = private)
        
        radiationMode;     % radiation modality
        identifier;       % upper case short notation of the current model in combination with the quantity used for optimization (e.g. LEM_RBExD) probably not needed
        bioOpt;           % boolean indicating biological optimization (=true) or physical optimization (=false)
        model;            % upper case short notation of the current model (e.g. LEM)
        quantityOpt;      % quantity used for optimizaiton
        quantityVis;      % quantity used per default for visualization
        description;      % short description of the biological model
        RBE;              % constant RBE
        
        % beamMixingModel  = 'ZaiderRossi';
    end
    
    % constant public properties which are visible outside of this class
    properties(Constant = true)
        AvailableModels                 = {'none','constRBE','MCN','WED','LEM','HEL'};   % cell array determines available models - if cell is deleted then the corersponding model can not be generated
        AvailableradiationModealities   = {'photons','protons','helium','carbon'};
        AvailableQuantitiesForOpt       = {'physicalDose','effect','RBExD'};
        
        AvailableAlphaXBetaX = {[0.036 0.024],    'prostate';
            [0.089 0.287],    'rectum and normal tissue';
            [0.55 0.05],      'head and neck MCN';
            [0.0499 0.0238],  'brain tissue';
            [0.1 0.05],       'default values';
            [0.1 0.005],      'default values'; %
            [0.0081 0.0033],  'LEM IV AB 2.45'; %
            [0.0030 0.0015],  'LEM IV AB 2'}; %
        
    end
    
    % constant private properties which are only visible within this class
    properties(Constant = true, Access = private)
        
        constRBE_protons = 1.1;
        constRBE_photons = 1;
        
        %McNamara variable RBE model for protons
        p0_MCN   = 0.999064;     % according to https://www.ncbi.nlm.nih.gov/pubmed/26459756
        p1_MCN   = 0.35605;
        p2_MCN   = 1.1012;
        p3_MCN   = -0.0038703;
        
        %Carabe Fernandez variable RBE model for protons
        p0_CAR   = 0.843; % http://www.tandfonline.com/doi/abs/10.1080/09553000601087176?journalCode=irab20
        p1_CAR   = 0.154;
        p2_CAR   = 2.686;
        p3_CAR   = 1.09;
        p4_CAR   = 0.006;
        p5_CAR   = 2.686;
        
        %Wedenberg variable RBE model for protons
        p0_WED   = 1; % https://www.ncbi.nlm.nih.gov/pubmed/22909391
        p1_WED   = 0.434;
        p2_WED   = 1;
        
        % Linear Scaling model for protons
        p_lamda_1_1          = 0.008; %0.008; % according to Malte Frese https://www.ncbi.nlm.nih.gov/pubmed/20382482 (FITTED for head and neck patients !)
        p_corrFacEntranceRBE = 0.5;   %[kev/mum]
        p_upperLETThreshold  = 30;    %[kev/mum]
        P_lowerLETThreshold  = 0.3;   %[kev/mum]
        
        % data driven parametrization of helium ions https://iopscience.iop.org/article/10.1088/0031-9155/61/2/888 
        p0_HEL = 1.36938e-1;   % values refer to the quadratic exponential fit f_QE
        p1_HEL = 9.73154e-3;
        p2_HEL = 1.51998e-2;
        
    end
    
    %% methods
    
    % private methods go here
    methods (Access = private)
                       
        function this = setBioModel(this)
            
            matRad_cfg = MatRad_Config.instance();
            
            this.RBE = NaN;
            % check for valid combinations 
            boolCHECK    = false;
            
            while ~boolCHECK
                
                switch this.radiationMode
                    
                    case {'photons'}
                        
                        setDefaultValues = false;
                        switch this.quantityOpt
                            
                            case {'physicalDose'}
                                if strcmp(this.model,'none')
                                    boolCHECK        = true;
                                    this.bioOpt      = false;
                                    this.quantityVis = 'physicalDose';
                                else
                                    setDefaultValues = true;
                                end
                                
                            case {'RBExD'}
                                if sum(strcmp(this.model,{'constRBE', 'none'})) == 1
                                    this.RBE         = this.constRBE_photons;
                                    boolCHECK        = true;
                                    this.bioOpt      = false;
                                    this.quantityVis = 'RBExD';
                                else
                                    setDefaultValues = true;
                                end
                                
                            case {'effect'}
                                if strcmp( this.model,'none')
                                    boolCHECK        = true;
                                    this.bioOpt      = true;
                                    this.quantityVis = 'RBExD';
                                else
                                    setDefaultValues = true;   
                                end
                                
                            otherwise
                                matRad_cfg.dispWarning('matRad: Invalid biological optimization quantity: %s; using physical dose instead. \n',this.quantityOpt);
                                this.quantityOpt = 'physicalDose';
                        end
                        
                        if setDefaultValues
                            matRad_cfg.dispWarning('matRad: Invalid biological model: %s; using "none" instead. \n',this.model);
                            this.model       = 'none';
                            this.quantityVis = 'physicalDose';
                            this.quantityOpt = 'physicalDose';
                        end
                        
                    case {'protons'}
                        
                        switch this.quantityOpt
                            
                            case {'physicalDose'}
                                if strcmp( this.model, 'none')
                                    boolCHECK           = true;
                                    this.bioOpt         = false;
                                    this.quantityVis    = 'physicalDose';
                                else
                                    matRad_cfg.dispWarning('matRad: Invalid biological model: %s; using "none" instead. \n',this.model);
                                    this.model  = 'none';
                                end
                                
                            case {'RBExD'}
                                if sum(strcmp( this.model, {'constRBE','MCN','WED'})) == 1
                                    boolCHECK           = true;
                                    this.bioOpt         = true;
                                    this.quantityVis    = 'RBExD';
                                    if sum(strcmp( this.model, {'constRBE'})) == 1
                                        this.RBE  = this.constRBE_protons;
                                        this.bioOpt     = false;
                                    end  
                                else
                                    matRad_cfg.dispWarning(['matRad: Invalid biological model: ' this.model  '; using constant RBE instead. \n']);
                                    this.model  = 'constRBE';
                                end
                                
                            case {'effect'}
                                if sum(strcmp( this.model, {'MCN','WED'})) == 1
                                    boolCHECK           = true;
                                    this.bioOpt         = true;
                                    this.quantityVis    = 'RBExD';
                                else
                                    matRad_cfg.dispWarning(['matRad: Invalid biological model: ' this.model  '; using MCN Model instead. \n']);
                                    this.model  = 'MCN';
                                end
                                
                            otherwise
                                matRad_cfg.dispWarning(['matRad: Invalid biological optimization quantity: ' this.quantityOpt  '; using "none" instead. \n']);
                                this.quantityOpt = 'physicalDose';
                        end
                       
                    case {'helium'}
                        switch this.quantityOpt
                            
                            case {'physicalDose'}
                                if strcmp( this.model, 'none')
                                    boolCHECK           = true;
                                    this.bioOpt         = false;
                                    this.quantityVis    = 'physicalDose';
                                else
                                    matRad_cfg.dispWarning(['matRad: Invalid biological model: ' this.model  '; using "none" instead. \n']);
                                    this.model  = 'none';
                                end
                                
                            case {'effect','RBExD'}
                                if sum(strcmp(this.model,{'LEM','HEL'})) > 0
                                    boolCHECK           = true;
                                    this.bioOpt         = true;
                                    this.quantityVis    = 'RBExD';
                                else
                                    matRad_cfg.dispWarning(['matRad: Invalid biological Model: ' this.model  '; using "none" instead. \n']);
                                    this.model = 'none';
                                end
                                
                            otherwise
                                matRad_cfg.dispWarning(['matRad: Invalid biological optimization quantity: ' this.quantityOpt  '; using "none" instead. \n']);
                                this.quantityOpt = 'physicalDose';
                        end
                        
                    case {'carbon'}
                        switch this.quantityOpt
                            
                            case {'physicalDose'}
                                if strcmp( this.model, 'none')
                                    boolCHECK           = true;
                                    this.bioOpt         = false;
                                    this.quantityVis    = 'physicalDose';
                                else
                                    matRad_cfg.dispWarning(['matRad: Invalid biological model: ' this.model  '; using "none" instead. \n']);
                                    this.model  = 'none';
                                end
                                
                            case {'effect','RBExD'}
                                if strcmp(this.model,'LEM')
                                    boolCHECK           = true;
                                    this.bioOpt         = true;
                                    this.quantityVis    = 'RBExD';
                                else
                                    matRad_cfg.dispWarning(['matRad: Invalid biological Model: ' this.model  '; using Local Effect Model instead. \n']);
                                    this.model = 'LEM';
                                end
                                
                            otherwise
                                matRad_cfg.dispWarning(['matRad: Invalid biological optimization quantity: ' this.quantityOpt  '; using "none" instead. \n']);
                                this.quantityOpt = 'physicalDose';
                        end
                                      
                    otherwise
                        matRad_cfg.dispWarning(['matRad: Invalid biological radiation mode: ' this.radiationMode  '; using photons instead. \n']);
                        this.radiationMode = 'photons';
                end
            end
           
            
            % check quantity for optimization
            if this.bioOpt 
                if sum(strcmp(this.quantityOpt,{'physicalDose','RBExD','effect'})) == 0
                    matRad_cfg.dispError(['matRad: Invalid quantity for optimization: ' this.quantityOpt  ]);
                end
            else
                if sum(strcmp(this.quantityOpt,{'physicalDose','RBExD'})) == 0
                    matRad_cfg.dispError(['matRad: Invalid quantity for optimization: ' this.quantityOpt  ]);
                end
            end
           
            
            % check quantity for visualization
            if this.bioOpt 
                if sum(strcmp(this.quantityVis,{'physicalDose','RBExD','effect'})) == 0
                    matRad_cfg.dispError(['matRad: Invalid quantity for visualization: ' this.quantityVis  ]);
                end
            else
                if sum(strcmp(this.quantityVis,{'physicalDose','RBExD'})) == 0
                    matRad_cfg.dispError(['matRad: Invalid quantity for visualization: ' this.quantityVis  ]);
                end
            end
            
            
            
        end % end of this = setBioModel(this)
    end % end off access methods private
    
    
    
    % public methods go here
    
    methods
        
        % default constructor
        function this = matRad_bioModel(sRadiationMode,sQuantityOpt, sModel)
            this.radiationMode = sRadiationMode;
            this.quantityOpt   = sQuantityOpt;       % setter checks for valid strings but not for valid combinations (e.g. photons_LEM
            this.model         = sModel;
            this               = setBioModel(this);
        end % end constructor
        
        
        % setter functions
        function this = set.radiationMode(this,value)
            if ischar(value) && sum(strcmp(value,{'photons','protons','helium','carbon'})) == 1
                this.radiationMode = value;
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('matRad: Cannot set radiation modality')
            end
        end
        
        
        function this = set.bioOpt(this,value)
            if islogical(value)
                this.bioOpt = value;
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('matRad: Cannot set bioOpt option')
            end
        end
        
        
        function this = set.model(this,value)
            if ischar(value)
                this.model = value;
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('matRad: Cannot set model option')
            end
        end
        
        function this = set.quantityOpt(this,value)
            if ischar(value)
                this.quantityOpt = value;
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('matRad: Cannot set quantityOpt option')
            end
        end
        
        function this = set.quantityVis(this,value)
            if ischar(value)
                this.quantityVis = value;
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('matRad: Cannot set quantityVis option')
            end
        end
        
        function this = set.description(this,value)
            if ischar(value)
                this.description = value;
            else
                matRad_cfg = MatRad_Config.instance();
                matRad_cfg.dispError('matRad: Cannot set description option')
            end
        end
        
        
        
        
        function [bixelAlpha,bixelBeta] = calcLQParameter(this,vRadDepths,baseDataEntry,mTissueClass,vAlpha_x,vBeta_x,vABratio)
            matRad_cfg = MatRad_Config.instance();
            
            bixelAlpha = NaN*ones(numel(vRadDepths),1);
            bixelBeta  = NaN*ones(numel(vRadDepths),1);
            
            % range shift
            depths = baseDataEntry.depths + baseDataEntry.offset;
            
            switch [this.radiationMode '_' this.model]
               
               case {'protons_constRBE'}
                  
               case {'protons_LSM'}
                  
                  bixelLET = matRad_interp1(depths,baseDataEntry.LET,vRadDepths);
                  bixelLET(isnan(bixelLET)) = 0;
                  
                  ix                 = this.p_lowerLETThreshold < bixelLET < this.p_upperLETThreshold;
                  
                  alpha_0            = vAlpha_x - (this.p_lamda_1_1 * this.p_corrFacEntranceRBE);
                  bixelAlpha(ixLSM)  = alpha_0(ix) + this.p_lamda_1_1 * bixelLET;
                  
                  if sum(ix) < length(bixelLET)
                     bixelAlpha(bixelLET > pln.bioParam.lowerLETThreshold) =  alpha_0(bixelLET > this.p_upperLETThreshold) + this.p_lamda_1_1 * this.p_upperLETThreshold;
                     bixelAlpha(bixelLET < pln.bioParam.lowerLETThreshold) =  alpha_0(bixelLET < this.p_lowerLETThreshold) + this.p_lamda_1_1 * this.p_lowerLETThreshold;
                  end
                  
                  bixelBeta        = vBeta_x;
                  
                  % TODO: assign normal tissue an RBE of 1.1
                  
               case {'protons_MCN'}
                  
                  
                  bixelLET = matRad_interp1(depths,baseDataEntry.LET,vRadDepths);
                  bixelLET(isnan(bixelLET)) = 0;
                  
                  RBEmax     = this.p0_MCN + ((this.p1_MCN * bixelLET )./ vABratio);
                  RBEmin     = this.p2_MCN + (this.p3_MCN  * sqrt(vABratio) .* bixelLET);
                  bixelAlpha = RBEmax    .* vAlpha_x;
                  bixelBeta  = RBEmin.^2 .* vBeta_x;
                  
               case {'protons_WED'}
                  
                  bixelLET = matRad_interp1(depths,baseDataEntry.LET,vRadDepths);
                  bixelLET(isnan(bixelLET)) = 0;
                  
                  RBEmax     = this.p0_WED + ((this.p1_WED * bixelLET )./ vABratio);
                  RBEmin     = this.p2_WED;
                  bixelAlpha = RBEmax    .* vAlpha_x;
                  bixelBeta  = RBEmin.^2 .* vBeta_x;
                  
               case {'helium_HEL'}
                  
                  % data-driven RBE parametrization of helium ions
                  % https://iopscience.iop.org/article/10.1088/0031-9155/61/2/888
                  
                  bixelLET = matRad_interp1(depths,baseDataEntry.LET,vRadDepths);
                  bixelLET(isnan(bixelLET)) = 0;
                  
                  % quadratic fit
                  %f_Q      = 8.53959e-4 .* bixelLET.^2;
                  %RBEmax_Q = 1 + 2.145e-1  + vABratio.^-1 .* f_Q;
                  % linear quadratic fit
                  %f_LQ      = 2.91783e-1*bixelLET - 9.525e-4*bixelLET.^2;
                  %RBEmax_LQ = 1 + ((1.42057e-1 + (vABratio.^-1)) .* f_LQ);
                  % linear exponential fit
                  %f_LE      = (2.965e-1 * bixelLET) .* exp(-4.90821e-3 * bixelLET);
                  %RBEmax_LE = 1 + ((1.5384e-1  + (vABratio.^-1)) .* f_LE);
                  
                  % quadratic exponential fit
                  f_QE      = (this.p1_HEL * bixelLET.^2) .* exp(-this.p2_HEL * bixelLET);
                  RBEmax_QE = 1 + ((this.p0_HEL  + (vABratio.^-1)) .* f_QE);

                  % the linear quadratic fit yielded the best fitting result
                  RBEmax = RBEmax_QE;
                  
                  RBEmin = 1; % no gain in using fitted parameters over a constant value of 1
                  
                  bixelAlpha = RBEmax    .* vAlpha_x;
                  bixelBeta  = RBEmin.^2 .* vBeta_x;
                  
               case {'carbon_LEM','helium_LEM'}
                  
                  numOfTissueClass = size(baseDataEntry(1).alpha,2);
                  
                  for i = 1:numOfTissueClass
                     bixelAlpha(mTissueClass==i) = matRad_interp1(depths,baseDataEntry.alpha(:,i),vRadDepths(mTissueClass==i));
                     bixelBeta(mTissueClass==i)  = matRad_interp1(depths,baseDataEntry.beta(:,i), vRadDepths(mTissueClass==i));
                  end
                  
               otherwise
                  
            end
            
        end
        
        

    end % end public methods
    
    
    methods(Static)
        
    end % end static public methods
    
    
end % end class definition