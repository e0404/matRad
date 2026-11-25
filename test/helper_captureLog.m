function [logEntries, varargout] = helper_captureLog(func)
% helper_captureLog runs func and returns the matRad log messages it produced
%
% The messages are recorded through MatRad_Config's keepLog switch, which
% stores them independent of the log level. Hence no console output is needed
% to check for warnings.
%
% call:
%   [logEntries, out1, out2, ...] = helper_captureLog(func)
%
% input:
%   func:       function handle without input arguments
%
% output:
%   logEntries: N x 2 cell array with the log type ('WARNING', 'ERROR', ...)
%               and the message of each log entry produced by func
%   varargout:  outputs of func

matRad_cfg = MatRad_Config.instance();
previousKeepLog = matRad_cfg.keepLog;
cleanup = onCleanup(@() helper_restoreKeepLog(matRad_cfg, previousKeepLog));

matRad_cfg.keepLog = true;
numPreviousEntries = size(matRad_cfg.messageLog, 1);

if nargout > 1
    [varargout{1:nargout - 1}] = func();
else
    func();
end

logEntries = matRad_cfg.messageLog(numPreviousEntries + 1:end, :);

end

function helper_restoreKeepLog(cfg, keepLog)
cfg.keepLog = keepLog;
end
