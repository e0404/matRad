function [ bioOptParam ] = matRad_getBioOptParam(sBioOptIdentifier,sRadiationMode)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% check for supported optimization methods using PHOTONS
if isequal(sRadiationMode,'photons')
    
    if isequal(sBioOptIdentifier,'physical')
        bioOptParam.Type             = 'physical';
    else
         warning(['matRad: Invalid biological optimization: ' sBioOptIdentifier ' for ' sRadiationMode '; using physical optimization instead'])
    end
    
% check for supported optimization methods using PROTONS    
elseif isequal(sRadiationMode,'protons')
    
    if isequal(sBioOptIdentifier,'physical')
        bioOptParam.Type             = 'physical';
        
    elseif isequal(sBioOptIdentifier,'const_RBExD') 
        bioOptParam.Type             = sBioOptIdentifier;
        bioOptParam.Model            = 'constantRBE';
        bioOptParam.Quantity         = 'RBExD';
        bioOptParam.constRBE         = 1.1;
        
    elseif isequal(sBioOptIdentifier,'LSM_effect')  
        bioOptParam.Type               = sBioOptIdentifier;
        bioOptParam.Model              = 'linear scaling model LSM';
        bioOptParam.Quantity           = 'effect';
        bioOptParam.lamda_1_1          = 0.008; %0.008; % according to Malte Frese https://www.ncbi.nlm.nih.gov/pubmed/20382482 (fitted for head and neck patients)
        bioOptParam.corrFacEntranceRBE = 0.5;   %[kev/mum]
        bioOptParam.upperLETThreshold  = 30;    %[kev/mum]
        bioOptParam.lowerLETThreshold  = 0.3;   %[kev/mum]
        
     elseif isequal(sBioOptIdentifier,'LSM_RBExD') 
        bioOptParam.Type               = sBioOptIdentifier;
        bioOptParam.Model              = 'linear scaling model LSM';
        bioOptParam.Quantity           = 'RBExD';
        bioOptParam.lamda_1_1          = 0.008; %0.008; % according to Malte Frese https://www.ncbi.nlm.nih.gov/pubmed/20382482 (fitted for head and neck patients)
        bioOptParam.corrFacEntranceRBE = 0.5;   %[kev/mum]
        bioOptParam.upperLETThreshold  = 30;    %[kev/mum]
        bioOptParam.lowerLETThreshold  = 0.3;   %[kev/mum]
        
    else
         warning(['matRad: Invalid biological optimization: ' sBioOptIdentifier ' for ' sRadiationMode '; using const_RBExD optimization instead'])
         % back up solution
        bioOptParam.Type             = 'const_RBExD';
        bioOptParam.Model            = 'constantRBE';
        bioOptParam.Quantity         = 'RBExD';
        bioOptParam.constRBE         = 1.1;
    end
    
% check for supported optimization methods using CARBON   
elseif isequal(sRadiationMode,'carbons')
    
    if isequal(sBioOptIdentifier,'physical')
        
        bioOptParam.Type             = 'physical';
        
    elseif isequal(sBioOptIdentifier,'LEMIV_effect')  
        
        bioOptParam.Type            = sBioOptIdentifier;
        bioOptParam.Model            = 'LEMIV';
        bioOptParam.Quantity         = 'effect';
        bioOptParam.BeamMixingModel  = 'ZaiderRossi';
        bioOptParam.dummy = 1;
        
    elseif isequal(sBioOptIdentifier,'LEMIV_RBExD') 
        
        bioOptParam.Type             = sBioOptIdentifier;
        bioOptParam.Model            = 'LEMIV';
        bioOptParam.Quantity         = 'RBExD';
        bioOptParam.BeamMixingModel  = 'ZaiderRossi';
        
    else
         warning(['matRad: Invalid biological optimization: ' sBioOptIdentifier ' for ' sRadiationMode '; using LEMIV_RBExD optimization instead'])
         % back up solution
        bioOptParam.Type             = sBioOptIdentifier;
        bioOptParam.Model            = 'LEMIV';
        bioOptParam.Quantity         = 'RBExD';
        bioOptParam.BeamMixingModel  = 'ZaiderRossi';
    end


end

