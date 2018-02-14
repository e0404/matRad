classdef matRad_bioModel
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  matRad_bioModel
    %  This class creates all required biological model parameters according to
    % a given radiation modatlity and a given bio model identifier string.
    %
    % constructor call
    %   pln.bioParam = matRad_bioModel(sRadiationMode,sIdentifier)
    %
    %   e.g. pln.bioParam = matRad_bioModel('protons','constRBE_RBExD')
    %
    % input
    %   sRadiationMode:     radiation modality 'photons' 'protons' 'carbon'
    %   sIdentifier:        string to denote a biological model along with the
    %                       quantity used for optimization
    %                       none_physicalDose: physical optimization;                              constRBE_RBExD; constant RBE of 1.1;
    %                       MCN_effect; McNamara-variable RBE model for protons (effect based)     MCN_RBExD; McNamara-variable RBE model for protons (RBExD) based
    %                       WED_effect; Wedenberg-variable RBE model for protons (effect based)    WED_RBExD; Wedenberg-variable RBE model for protons (RBExD) based
    %                       LEM_effect: effect-based optimization;                                 LEM_RBExD: optimization of RBE-weighted dose
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
        AvailableModels                 = {'none','constRBE','MCN','WED','LEM'};   % cell array determines available models - if cell is deleted then the corersponding model can not be generated
        AvailableradiationModealities   = {'photons','protons','helium','carbon'};
        AvailableQuantitiesForOpt       = {'physicalDose','effect','RBExD'};
        
        AvailableAlphaXBetaX       = {[0.036 0.024],    'prostate';
            [0.089 0.287],    'rectum and normal tissue';
            [0.55 0.05],      'head and neck MCN';
            [0.0499 0.0238],  'brain tissue';
            [0.1 0.05],       'default values';
            [0.1 0.005],      'default values'}; %
        
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
        
    end
    
    %% methods
    
    % private methods go here
    methods (Access = private)
                       
        function this = setBioModel(this)
            
            this.RBE = NaN;
            % check for valid combinations 
            boolCHECK    = false;
            
            while ~boolCHECK
                
                switch this.radiationMode
                    case {'photons'}
                        
                        switch this.quantityOpt
                            
                            case {'physicalDose'}
                                if strcmp( this.model, 'none')
                                    boolCHECK           = true;
                                    this.bioOpt         = false;
                                    this.quantityVis    = 'physicalDose';
                                else
                                    matRad_dispToConsole(['matRad: Invalid biological model: ' this.model  '; using "none" instead. \n'],[],'warning');
                                    this.model  = 'none';
                                end
                                
                            case {'RBExD'}
                                if sum(strcmp(this.model, {'constRBE', 'none'})) == 1
                                    this.RBE            = this.constRBE_photons;
                                    boolCHECK           = true;
                                    this.bioOpt         = false;
                                    this.quantityVis    = 'RBExD';
                                else
                                    matRad_dispToConsole(['matRad: Invalid biological model: ' this.model  '; using constant RBE instead. \n'],[],'warning');
                                    this.model = 'constRBE';
                                end
                                
                            case {'effect'}
                                if strcmp( this.model, 'none')
                                    boolCHECK           = true;
                                    this.bioOpt         = true;
                                    this.quantityVis    = 'RBExD';
                                else
                                    matRad_dispToConsole(['matRad: Invalid biological model: ' this.model  '; using "none" instead. \n'],[],'warning');
                                    this.model  = 'none';
                                end
                                
                            otherwise
                                matRad_dispToConsole(['matRad: Invalid biological optimization quantity: ' this.quantityOpt  '; using physical dose instead. \n'],[],'warning');
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
                                    matRad_dispToConsole(['matRad: Invalid biological model: ' this.model  '; using "none" instead. \n'],[],'warning');
                                    this.model  = 'none';
                                end
                                
                            case {'RBExD'}
                                if sum(strcmp( this.model, {'constRBE','MCN','WED'})) == 1
                                    boolCHECK           = true;
                                    this.bioOpt         = false;
                                    this.quantityVis    = 'RBExD';
                                    if sum(strcmp( this.model, {'constRBE'})) == 1
                                        this.RBE  = this.constRBE_protons;
                                    end  
                                else
                                    matRad_dispToConsole(['matRad: Invalid biological model: ' this.model  '; using constant RBE instead. \n'],[],'warning');
                                    this.model  = 'constRBE';
                                end
                                
                            case {'effect'}
                                if sum(strcmp( this.model, {'MCN','WED'})) == 1
                                    boolCHECK           = true;
                                    this.bioOpt         = true;
                                    this.quantityVis    = 'RBExD';
                                else
                                    matRad_dispToConsole(['matRad: Invalid biological model: ' this.model  '; using MCN Model instead. \n'],[],'warning');
                                    this.model  = 'MCN';
                                end
                                
                            otherwise
                                matRad_dispToConsole(['matRad: Invalid biological optimization quantity: ' this.quantityOpt  '; using "none" instead. \n'],[],'warning');
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
                                    matRad_dispToConsole(['matRad: Invalid biological model: ' this.model  '; using "none" instead. \n'],[],'warning');
                                    this.model  = 'none';
                                end
                            otherwise
                                matRad_dispToConsole(['matRad: Invalid biological optimization quantity: ' this.quantityOpt  '; using "none" instead. \n'],[],'warning');
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
                                    matRad_dispToConsole(['matRad: Invalid biological model: ' this.model  '; using "none" instead. \n'],[],'warning');
                                    this.model  = 'none';
                                end
                                
                            case {'effect','RBExD'}
                                if strcmp(this.model,'LEM')
                                    boolCHECK           = true;
                                    this.bioOpt         = true;
                                    this.quantityVis    = 'RBExD';
                                else
                                    matRad_dispToConsole(['matRad: Invalid biological Model: ' this.model  '; using Local Effect Model instead. \n'],[],'warning');
                                    this.model = 'LEM';
                                end
                                
                            otherwise
                                matRad_dispToConsole(['matRad: Invalid biological optimization quantity: ' this.quantityOpt  '; using "none" instead. \n'],[],'warning');
                                this.quantityOpt = 'physicalDose';
                        end
                        
                    case {'helium'}
                     
                        warning('not yet implemented');  
                        
                    otherwise
                        matRad_dispToConsole(['matRad: Invalid biological radiation mode: ' this.radiationMode  '; using photons instead. \n'],[],'warning');
                        this.radiationMode = 'photons';
                end
            end
           
            
            % check quantity for optimization
            if this.bioOpt 
                if sum(strcmp(this.quantityOpt,{'physicalDose','RBExD','effect'})) == 0
                    matRad_dispToConsole(['matRad: Invalid quantity for optimization: ' this.quantityOpt  ],[],'error');
                end
            else
                if sum(strcmp(this.quantityOpt,{'physicalDose','RBExD'})) == 0
                    matRad_dispToConsole(['matRad: Invalid quantity for optimization: ' this.quantityOpt  ],[],'error');
                end
            end
           
            
            % check quantity for visualization
            if this.bioOpt 
                if sum(strcmp(this.quantityVis,{'physicalDose','RBExD','effect'})) == 0
                    matRad_dispToConsole(['matRad: Invalid quantity for visualization: ' this.quantityVis  ],[],'error');
                end
            else
                if sum(strcmp(this.quantityVis,{'physicalDose','RBExD'})) == 0
                    matRad_dispToConsole(['matRad: Invalid quantity for visualization: ' this.quantityVis  ],[],'error');
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
                error('matRad: Cannot set radiation modality')
            end
        end
        
        
        function this = set.bioOpt(this,value)
            if islogical(value)
                this.bioOpt = value;
            else
                error('matRad: Cannot set bioOpt option')
            end
        end
        
        
        function this = set.model(this,value)
            if ischar(value)
                this.model = value;
            else
                error('matRad: Cannot set model option')
            end
        end
        
        function this = set.quantityOpt(this,value)
            if ischar(value)
                this.quantityOpt = value;
            else
                error('matRad: Cannot set quantityOpt option')
            end
        end
        
        function this = set.quantityVis(this,value)
            if ischar(value)
                this.quantityVis = value;
            else
                error('matRad: Cannot set quantityVis option')
            end
        end
        
        function this = set.description(this,value)
            if ischar(value)
                this.description = value;
            else
                error('matRad: Cannot set description option')
            end
        end
        
        
        
        
        function [bixelAlpha,bixelBeta] = calcLQParameter(this,vRadDepths,baseDataEntry,mTissueClass,vAlpha_x,vBeta_x,vABratio)
            
            
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
                    
                case {'carbon_LEM'}
                    
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