function statusAll = runMatRadTests(varargin)
%% This file runs the complete matRad test suite.


CI_MODE = strcmpi(getenv('CONTINUOUS_INTEGRATION'),'true') || strcmp(getenv('CI'),'true');
isJenkins = ~isempty(getenv('JENKINS_URL'));

%% Set path
addpath(fullfile(pwd,'..','src'));
addpath(fullfile(pwd,'suites'));

%% Select functions to runs
suite = @ACID;
allTests = 1:numel(suite(0));

%% Prepare environment
if strcmpi(getEnvironment(), 'Octave')
  % Ensure that paging is disabled
  % https://www.gnu.org/software/octave/doc/interpreter/Paging-Screen-Output.html
  more off
end

%% Run tests
status = testHeadless('testFunctionIndices', allTests,...
                     'testsuite',           suite, varargin{:});

if isJenkins
    makeTapReport(status, 'stream', 'results.test.tap');
    makeTravisReport(status, 'stream', 'results.test.md');
end

nErrors = makeTravisReport(status);

%% Calculate exit code
if CI_MODE
    exit(nErrors);
end
