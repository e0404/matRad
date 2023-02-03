function dij = matRad_calcCombiDose(ct,stf,pln,cst,CalcDoseDirect)
% data structure handling and dij calculation for Joint optimization and
% spatiotemporal plans 
% 
% call
%   dij = matRad_calcCombiDose(ct,stf,pln,cst,CalcDoseDirect)
%
% input
%   ct :            matRad ct struct
%   stf:            matRad stf stuct
%   cst:            matRad cst struct
%   pln:            matRad pln struct
%   CalcDoseDirect: (optional) boolian to bypass dose influence matrix
%                   computation and directly calculate dose; only makes
%                   sense in combination with matRad_calcDoseDirect.m
%
% output
%   dij:            matRad dij struct
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    matRad_cfg = MatRad_Config.instance();

    if (strcmp(pln.radiationMode, 'MixMod'))

        matRad_cfg.dispInfo(' PROJECT MIXED MOD HAS BEEN ACTIVATED. AMDG. !! \n \n');
   
         stfModalities = {stf(:).radiationMode};
         radiationModalities = {pln.originalPlans.radiationMode};
         dij_fields = [];
         dijt = [];
         
         for k = 1:pln.numOfModalities

            currStf  = stf(strcmp(stfModalities,radiationModalities{k}));
            currPln = pln.originalPlans(k);
            currPln.bioParam = pln.bioParam(k);
            

            if strcmp(radiationModalities{k},'photons')
                dijt = [dijt, {matRad_calcPhotonDose(ct,currStf,currPln,cst,CalcDoseDirect)}];
            elseif strcmp(radiationModalities{k},'protons') || strcmp(radiationModalities{k},'carbon')|| strcmp(radiationModalities{k},'helium')
                dijt = [dijt, {matRad_calcParticleDose(ct,currStf,currPln,cst,CalcDoseDirect)}];
            end
   
             dij_fields = [dij_fields; fieldnames(dijt{k})];
         end

         [~, idx] = unique(dij_fields);
         idx = setdiff(1:numel(dij_fields),idx);

         commonFields = dij_fields(idx);
         commonFields(:,2) = cell(size(commonFields));
         commonFields = commonFields';

         dij = struct(commonFields{:});

         dij_fieldnames = fieldnames(dij);

         for k=1:size(dij_fieldnames,1)
            currField = dij_fieldnames{k};
            dijt_fields = [cellfun(@(x) x.(currField), dijt, 'UniformOutput', false)];
            if isequaln(dijt_fields{:})
                dij.(currField) = dijt_fields{1};
            end
         end

        fieldsName = 'numOfBeams';
        if isfield(dij, fieldsName)
            fieldValues = [cellfun(@(x) x.(fieldsName), dijt, 'UniformOutput', false)];
            dij.(fieldsName) = [fieldValues{:}];
        end

        fieldsName = 'totalNumOfBixels';
        if isfield(dij, fieldsName)
            fieldValues = [cellfun(@(x) x.(fieldsName), dijt, 'UniformOutput', false)];
            dij.(fieldsName) = sum([fieldValues{:}]);
        end

        fieldsName = 'totalNumOfRays';
        if isfield(dij, fieldsName)
            fieldValues = [cellfun(@(x) x.(fieldsName), dijt, 'UniformOutput', false)];
            dij.(fieldsName) = sum([fieldValues{:}]);
        end

        fieldsName = 'numOfRaysPerBeam';
        if isfield(dij, fieldsName)
            fieldValues = [cellfun(@(x) x.(fieldsName), dijt, 'UniformOutput', false)];
            dij.(fieldsName) = [fieldValues{:}];
        end

        spareStruct = 0;
        for i=1:size(cst,1)
            if ~isempty(cst{i,6}{1}) && numel(unique([cst{i,6}{1}.alphaX],'stable'))>1
                spareStruct = i;
            end
        end
        if spareStruct >0 
            dij.spareStruct = spareStruct;
            dij.alphaCubes = 2;
        else 
                dij.alphaCubes = 1;
        end


        for k=1:pln.numOfModalities
            if ~isfield (pln.propOpt(k),'spatioTemp')
                pln.propOpt(k).spatioTemp = 0;
                pln.propOpt(k).STScenarios = 1;
                pln.propOpt(k).STfractions = pln.originalPlans(k).numOfFractions;
            else
                if ~isfield(pln.propOpt(k), 'STScenarios') || pln.propOpt(k).spatioTemp == 0
                    pln.propOpt(k).STScenarios = 1;
                end
                if ~isfield(pln.propOpt(k),'STfractions')                           % adjust not fully divisible number of fraction .... later
                    pln.propOpt(k).STfractions = floor(pln.originalPlans(k).numOfFractions/pln.propOpt(k).STScenarios)* ones(pln.propOpt(k).STScenarios,1);
                end
            end
        end
        
        dij.spatioTemp = [pln.propOpt.spatioTemp];
        dij.numOfSTScen = [pln.propOpt.STScenarios];
        dij.numOfSTFractions = [pln.propOpt.STfractions];
        
        dij.numOfModalities = pln.numOfModalities;
        dij.totalNumOfFractions = pln.numOfFractions;

        dij.dualIrradiation = false;

        dij.original_Dijs = [dijt];


         matRad_cfg.dispInfo('SUCCESS. I Dij It  !! \n');


    else
    
        if (strcmp(pln.radiationMode, 'photons'))
            dij = matRad_calcPhotonDose(ct,stf,pln,cst,CalcDoseDirect);
        elseif (strcmp(pln.radiationMode, 'protons') || strcmp(pln.radiationMode, 'carbon') || strcmp(pln.radiationMode, 'helium'))
            dij = matRad_calcParticleDose(ct,stf,pln,cst,CalcDoseDirect);
        end
    end


end