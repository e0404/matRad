function found = helper_logContains(logEntries, logType, pattern)
% helper_logContains checks whether a log entry of the given type contains
% the given text
%
% call:
%   found = helper_logContains(logEntries, logType, pattern)
%
% input:
%   logEntries: N x 2 cell array of log entries, see helper_captureLog
%   logType:    log type, e.g. 'WARNING' or 'ERROR'
%   pattern:    text that has to be part of the message
%
% output:
%   found:      true if any log entry of type logType contains pattern

if isempty(logEntries)
    found = false;
    return
end

messages = logEntries(strcmpi(logEntries(:, 1), logType), 2);
found = any(~cellfun(@isempty, strfind(messages, pattern)));

end
