function [machine] = matRad_getAlphaBetaCurves(machine,varargin)
% matRad alpha beta curve calculation tool
%
% call
%   machine = matRad_getAlphaBetaCurves(machine)
%   machine = matRad_getAlphaBetaCurves(machine,cst,modelName,overrideAB)
% Example full call for protons
%   machine = matRad_getAlphaBetaCurves(machine,pln,cst,'MCN','override')
% input
%   machine:                matRad machine file to change
%   varargin (optional):    cst:        matRad cst struct (for custom alpha/beta,
%                                       otherwise default is alpha=0.1, beta=0.05;
%                           modelName:  specify RBE modelName
%                   	    overrideAB: calculate new alpha beta even if available
%                                       and override
%
% output
%   machine:                updated machine file with alpha/beta curves
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2021 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg =  MatRad_Config.instance();

overrideAB = false;
if ~isempty(varargin)
    for i = 1:nargin-1
        if iscell(varargin{i})
            cst = varargin{i};
        elseif ischar(varargin{i})
            if contains(varargin{i},'override')
                if  isfield(machine.data,'alphaX')
                    overrideAB = true;
                end
            else
                modelName = varargin{i};
            end
        end
    end
end

if ~contains(machine.meta.radiationMode,{'protons','helium'}) || (exist('modelName','var') && contains(modelName,'LEM'))
    matRad_cfg.dispWarning('alphaBetaCurves cannot be calculated for LEM model, skipping...');
elseif ~isfield(machine.data,'alphaX') || overrideAB
    %% save original fields with suffix if available
    if overrideAB
        fieldnames = {'alphaX','betaX','alphaBetaRatio','alpha','beta'};
        for k=1:numel(fieldnames)
            [machine.data.([fieldnames{k} '_org'])] = deal(machine.data.(fieldnames{k}));
        end
        machine.data = rmfield(machine.data,fieldnames);
    end
    
    %% Set biological models
    if ~exist('modelName','var')
        if contains(machine.meta.radiationMode,'proton')
            modelName = 'MCN'; % Use McNamara model as default for protons
        elseif contains(machine.meta.radiationMode,'helium')
            modelName = 'HEL';
        end
    end
    pln.bioParam = matRad_bioModel(machine.meta.radiationMode,'RBExD',modelName);
    
    %% get unique combintions of alpha/beta from cst or use default alpha/beta values
    if ~exist('cst','var')
        ab(1,1) = 0.1;
        ab(1,2) = 0.05;
    else
        ab = cellfun(@(f)cst{f,5}.alphaX,num2cell(1:10)');
        ab(:,2) = cellfun(@(f)cst{f,5}.betaX,num2cell(1:10)');
    end
    ab = unique(ab,'rows');
    
    % save alpha/beta values for each energy in machine.data
    [machine.data(:).alphaX] = deal(ab(:,1));
    [machine.data(:).betaX] = deal(ab(:,2));
    [machine.data(:).alphaBetaRatio] = deal(ab(:,1) ./ ab(:,2));
    
    % calculate alpha/beta curves for each energy in machine.data
    for j = 1:length(machine.data)
        depths = machine.data(j).depths + machine.data(j).offset;
        voxelsOnes = ones(numel(depths),1);
        
        machine.data(j).alpha = zeros(numel(voxelsOnes),size(ab,1));
        machine.data(j).beta = zeros(numel(voxelsOnes),size(ab,1));
        
        for k = 1:length(machine.data(j).alphaX)
            [machine.data(j).alpha(:,k),machine.data(j).beta(:,k)] = pln.bioParam.calcLQParameter(depths,machine.data(j),voxelsOnes,...
                machine.data(j).alphaX(k)*voxelsOnes,machine.data(j).betaX(k)*voxelsOnes,machine.data(j).alphaBetaRatio(k)*voxelsOnes);
        end
    end
    
else
    matRad_cfg.dispWarning('Basedata already contains alpha/beta curves. Use "overrideAB"-flag to override.')
end

end