function cleanup = helper_setStrictSelection(isStrict)
% helper_setStrictSelection sets MatRad_Config's strictSelection switch and
% restores the previous value when the returned object is cleared
%
% call:
%   cleanup = helper_setStrictSelection(isStrict)
%
% input:
%   isStrict:   value for MatRad_Config.instance().strictSelection
%
% output:
%   cleanup:    onCleanup object, keep it until the end of the test

matRad_cfg = MatRad_Config.instance();
previousValue = matRad_cfg.strictSelection;
cleanup = onCleanup(@() helper_restoreStrictSelection(matRad_cfg, previousValue));
matRad_cfg.strictSelection = isStrict;

end

function helper_restoreStrictSelection(cfg, value)
cfg.strictSelection = value;
end
