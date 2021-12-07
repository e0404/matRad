function [machine] = matRad_getAlphaBetaCurves(machine,varargin)
% matRad alpha beta curve calculation tool
%
% call
%   machine.data = matRad_getAlphaBetaCurves(machine.data,pln,cst)
%
% input
%   machine.data:       matRad machine file to change
%   varargin:       (optional) matRad cst struct and/or specify RBE modelName
%                   to use for computation
%
% output
%   machine.data:       updated machine file with alpha/beta curves
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

if contains(machine.meta.radiationMode,{'car','hel'})
   matRad_cfg.dispError('alphaBetaCurves cannot be calculated for carbon or helium ions.'); 
end

if ~isempty(varargin)
    for i = 1:nargin-1
        if iscell(varargin{i})
            cst = varargin{i};
        elseif isstr(varargin{i})
            modelName = varargin{i};
        end
    end
end

%% Set biological models
pln.radiationMode   = machine.meta.radiationMode;     % either photons / protons / carbon
if ~exist('modelName','var')
    if contains(pln.radiationMode,'proton')
        modelName = 'MCN'; % Use McNamara model as default for protons
    end
end
pln.bioParam = matRad_bioModel(pln.radiationMode,'RBExD',modelName);

% save original fields with suffix if available
if isfield(machine.data,'alphaX')
    fieldnames = {'alphaX','betaX','alphaBetaRatio','alpha','beta'};
    for k=1:numel(fieldnames)
        [machine.data.([fieldnames{k} '_org'])] = deal(machine.data.(fieldnames{k}));
    end
    machine.data = rmfield(machine.data,fieldnames);
end

%% get unique combinations of alpha/beta from cst or use default alpha/beta values
if ~exist('cst','var')
    ab(1,1) = 0.1;
    ab(1,2) = 0.05;
else
    for i = 1:size(cst,1)
        ab(i,1) = cst{i,5}.alphaX;
        ab(i,2) = cst{i,5}.betaX;
    end
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

end