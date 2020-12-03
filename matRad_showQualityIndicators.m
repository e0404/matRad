function matRad_showQualityIndicators(qi)
% matRad display of quality indicators as table
% 
% call
%   matRad_showQualityIndicators(qi)
%
% input
%   qi: result struct from matRad_calcQualityIndicators
%
% output
%   graphical display of quality indicators in table form   
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team. 
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

[env, vStr] = matRad_getEnvironment();
    
% Create the column and row names in cell arrays 
rnames = {qi.name};
qi = rmfield(qi,'name');
cnames = fieldnames(qi);
for i = 1:numel(cnames)
    ix = find(cnames{i}(4:end) == '_');
    if ~isempty(ix)
        cnames{i}(ix+3) = '.';
    end
end

%To avoid parse error in octave, replace empty qi values with '-'
qi = (squeeze(struct2cell(qi)))';
qiEmpty = cellfun(@isempty,qi);
qi(qiEmpty) = {'-'};

%since uitable is only available in newer octave versions, we try and catch
try
    % Create the uitable
    table = uitable(gcf,'Data',qi,...
        'ColumnName',cnames,...
        'RowName',rnames,'ColumnWidth',{70});
    
    % Layout
    pos = get(gca,'position');
    set(table,'units','normalized','position',pos)
    axis off
catch ME
    matRad_cfg.dispWarning('The uitable function is not implemented in %s v%s.',env,vStr);
end
