function varargout = UnitTest(varargin)
% UnitTest executes the unit tests for this application, and can be called 
% either independently (when testing just the latest version) or via 
% UnitTestHarness (when testing for regressions between versions).  Either 
% two or three input arguments can be passed to UnitTest as described 
% below.
%
% The following variables are required for proper execution: 
%   varargin{1}: string containing the path to the main function. This is
%       only required if UnitTest is to be operated across multiple file 
%       names; otherwise it can be an empty string.
%   varargin{2}: variable containing the test data.  this can be a string
%       to a file, a structure, or any MATLAB variable.
%   varargin{3} (optional): variable containing reference data to be used
%       for comparison.  If not provided, it is assumed that this version
%       is the reference and therefore all comparison tests will "Pass".
%
% The following variables are returned upon succesful completion when input 
% arguments are provided. If no return variables are given, the results
% will be printed to stdout.
%   varargout{1}: cell array of strings containing preamble text that
%       summarizes the test, where each cell is a line. This text will
%       precede the results table in the report.
%   varargout{2}: n x 3 cell array of strings containing the test ID in
%       the first column, name in the second, and result (Pass/Fail or 
%       numerical values typically) of the test in the third.
%   varargout{3}: cell array of strings containing footnotes referenced by
%       the tests, where each cell is a line.  This text will follow the
%       results table in the report.
%   varargout{4} (optional): variable containing reference data created by 
%       executing this version.  This variable can be passed back into 
%       subsequent executions of UnitTest as varargin{3} to compare results
%       between versions (or to a priori validated reference data).
%
% The following variables are returned when no input arguments are
% provided (required only if called by UnitTestHarness):
%   varargout{1}: string containing the application name (with .m 
%       extension)
%   varargout{2}: string containing the path to the version application 
%       whose results will be used as reference
%   varargout{3}: 1 x n cell array of strings containing paths to the other 
%       applications which will be tested
%   varargout{4}: 2 x m cell array of strings containing the name of each 
%       test suite (first column) and path to the test data (second column)
%   varargout{5}: string containing the path and name of report file (will 
%       be appended by _R201XX.md based on the MATLAB version)
%
% Below is an example of how this function is used:
%
%   % Start profiler
%   profile on -history
%
%   % Execute unit test, printing the test results to stdout
%   UnitTest('', load('../test_data/InputData.mat'), ...
%       load('../test_data/OutputData.mat'));
%
%   % Stop and view the profiler results
%   profile viewer
%
%   % Execute unit test, storing the test results
%   [preamble, table, footnotes] = UnitTest('', ...
%       load('../test_data/InputData.mat'), ...
%       load('../test_data/OutputData.mat'));
%
%   % Execute unit test again but without reference data, this time storing 
%   % the output from UnitTest as a new reference file
%   [preamble, table, footnotes, newreference] = ...
%       UnitTest('', load('../test_data/InputData.mat'));
%
% Author: Mark Geurts, mark.w.geurts@gmail.com
% Copyright (C) 2015 University of Wisconsin Board of Regents
%
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the  
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along 
% with this program. If not, see http://www.gnu.org/licenses/.

%% Return Application Information
% If UnitTest was executed without input arguments
if nargin == 0
    
    % Declare the application filename
    varargout{1} = '';

    % Declare current version directory
    varargout{2} = '';

    % Declare prior version directories
    varargout{3} = {};

    % Declare location of test data. Column 1 is the name of the 
    % test suite, column 2 is the absolute path to the file(s)
    varargout{4} = {};

    % Declare name of report file (will be appended by _R201XX.md based on 
    % the MATLAB version)
    varargout{5} = '';
    
    % Return to invoking function
    return;
end

%% Initialize Unit Testing
% Initialize static test result text variables
pass = 'Pass';
fail = 'Fail';
unk = 'N/A'; %#ok<NASGU>

% Initialize preamble text
preamble = {
    '| Input Data | Value |'
    '|------------|-------|'
};

% Initialize results cell array
results = cell(0,3);

% Initialize footnotes cell array
footnotes = cell(0,1);

% Initialize reference structure
if nargout == 4
    reference = struct;
end

%% TEST 1/2/3/4: Analytic 1D Gamma
%
% DESCRIPTION: This unit test creates a large analytic reference dataset, 
%   modifies it by a known amount, and executes CalcGamma.  The gamma
%   result is then compared to a known solution
%
% RELEVANT REQUIREMENTS: none
%
% INPUT DATA: varargin{3}.gammaA
%
% CONDITION A (+): When computed using identical criteria, CalcGamma
%   produces an identical array to the reference with GPU Calculation
%
% CONDITION B (+): When computed using identical criteria, CalcGamma
%   produces an identical array to the reference with CPU Calculation
%
% CONDITION C (-): When computed using modified criteria, CalcGamma
%   produces a different result

% Add path to Event() function (to test Event calls function)
addpath('../test_data/');

% Declare gamma criteria
percent = 3;
dta = 0.3;
local = 1;

% Declare size of data arrays
n = 2000;

% Create reference structure (sine function)
ref.start = 0;
ref.width = pi / 1000;
ref.data = sin(ref.start:ref.width:...
    ref.start + ref.width * n) + 2;

% Declare modification values (to be applied to target data)
shift = 0.1;
scale = 1.01;

% Modify reference data, saving to target structure
target.start = ref.start + shift;
target.width = ref.width;
target.data = ref.data * scale;

% Execute in try/catch statement
try

    % Execute CalcGamma using test dataset
    t = tic;
    gamma = CalcGamma(ref, target, percent, dta, 'local', local);
    gtime = sprintf('%0.1f sec', toc(t));
    gpf = pass;
    cpf = pass;

    % If reference data exists
    if nargin == 3
        
        % If current value does not equal the reference
        if max(abs(gamma - varargin{3}.gammaA)) > 1e-5
            gpf = fail;
        end
        
        % Execute CalcGamma again using CPU
        t = tic;
        gamma = CalcGamma(ref, target, percent, dta, 'cpu', 1, ...
            'local', local);
        ctime = sprintf('%0.1f sec', toc(t));
        
        % If CPU value does not equal the reference
        if max(abs(gamma - varargin{3}.gammaA)) > 1e-5
            cpf = fail;
        end
        
        % Execute CalcGamma again, this time using modified criteria
        gamma = CalcGamma(ref, target, percent, dta * 2, ...
            'local', local);
        
        % If modified value equals the reference, the test fails
        if max(abs(gamma - varargin{3}.gammaA)) < 1e-5
            gpf = fail;
        end
        
        % Execute CalcGamma again, this time using modified criteria
        gamma = CalcGamma(ref, target, percent * 2, dta, ...
            'local', local);
        
        % If modified value equals the reference, the test fails
        if max(abs(gamma - varargin{3}.gammaA)) < 1e-5
            gpf = fail;
        end
        
    % Otherwise, set reference    
    elseif nargout == 4
        reference.gammaA = gamma;
        gpf = pass;
        cpf = pass;
    end
    
catch
    gpf = fail;
    cpf = fail;
end

% Add the dataset size to the preamble
preamble{length(preamble)+1} = ...
    sprintf('| Analytic 1D Gamma Array | %i x %i  |', size(gamma));

% Add header to results table
results{size(results,1)+1,1} = 'ID';
results{size(results,1),2} = 'Test Case';
results{size(results,1),3} = 'Result';

% Add GPU computation result
results{size(results,1)+1,1} = '1';
results{size(results,1),2} = 'Analytic 1D Gamma GPU Result within 1e-5';
results{size(results,1),3} = gpf;

% Add computation time
results{size(results,1)+1,1} = '2';
results{size(results,1),2} = 'Analytic 1D Gamma GPU Computation Time';
results{size(results,1),3} = gtime;

% Add CPU computation result
results{size(results,1)+1,1} = '3';
results{size(results,1),2} = 'Analytic 1D Gamma CPU Result within 1e-5';
results{size(results,1),3} = cpf;

% Add computation time
results{size(results,1)+1,1} = '4';
results{size(results,1),2} = 'Analytic 1D Gamma CPU Computation Time';
results{size(results,1),3} = ctime;

%% TEST 5/6/7/8: Beam Profile 1D Gamma
%
% DESCRIPTION: This unit test creates a real world profile dataset 
%   and executes CalcGamma.  The gamma result is then compared to a known 
%   solution
%
% RELEVANT REQUIREMENTS: none
%
% INPUT DATA: varargin{2}.referenceB, varargin{2}.targetB, 
%   varargin{3}.gammaB
%
% CONDITION A (+): When computed using identical criteria, CalcGamma
%   produces an identical array to the reference with GPU Calculation
%
% CONDITION B (+): When computed using identical criteria, CalcGamma
%   produces an identical array to the reference with CPU Calculation
%
% CONDITION C (-): When computed using modified criteria, CalcGamma
%   produces a different result

% Declare gamma criteria
percent = 2;
dta = 0.1;
local = 0;

% Load reference and target data from varargin
ref = varargin{2}.referenceB;
target = varargin{2}.targetB;

% Determine location and value of maximum
[C, I] = max(target.data);
    
% Search left side for half-maximum value
for x = 1:I-1
    if target.data(x) == C/2
        left = target.start + target.width * (x-1);
        break;
    elseif target.data(x) < C/2 && target.data(x+1) > C/2
        left = interp1(target.data(x:x+1), target.start + ...
            target.width * (x-1:x), C/2);
        break;
    end
end

% Search right side for half-maximum value
for x = I:size(target.data,1)-1
    if target.data(x) == C/2
        right = target.start + target.width * (x-1);
        break;
    elseif target.data(x) > C/2 && target.data(x+1) < C/2
        right = interp1(target.data(x:x+1), target.start + ...
            target.width * (x-1:x), C/2);
        break;
    end
end     

% Center target profile on FWHM
target.start = target.start - (right + left)/2;

% Execute in try/catch statement
try

    % Execute CalcGamma using test dataset
    t = tic;
    gamma = CalcGamma(ref, target, percent, dta, 'local', local);
    gtime = sprintf('%0.1f sec', toc(t));
    gpf = pass;
    cpf = pass;

    % If reference data exists
    if nargin == 3
        
        % If current value does not equal the reference
        if max(abs(gamma - varargin{3}.gammaB)) > 1e-5
            gpf = fail;
        end
        
        % Execute CalcGamma again using CPU
        t = tic;
        gamma = CalcGamma(ref, target, percent, dta, 'cpu', 1, ...
            'local', local);
        ctime = sprintf('%0.1f sec', toc(t));
        
        % If CPU value does not equal the reference
        if max(abs(gamma - varargin{3}.gammaB)) > 1e-5
            cpf = fail;
        end
        
        % Execute CalcGamma again, this time using modified criteria
        gamma = CalcGamma(ref, target, percent, dta, ...
            'local', local, 'limit', 0.1);
        
        % If modified value equals the reference, the test fails
        if max(abs(gamma - varargin{3}.gammaB)) < 1e-5
            gpf = fail;
        end
        
        % Execute CalcGamma again, this time using modified criteria
        gamma = CalcGamma(ref, target, percent, dta, ...
            'local', local, 'refval', 10);
        
        % If modified value equals the reference, the test fails
        if max(abs(gamma - varargin{3}.gammaB)) < 1e-5
            gpf = fail;
        end
        
    % Otherwise, set reference    
    elseif nargout == 4
        reference.gammaB = gamma;
        gpf = pass;
        cpf = pass;
    end
    
catch
    gpf = fail;
    cpf = fail;
end

% Add the dataset size to the preamble
preamble{length(preamble)+1} = ...
    sprintf('| Beam Profile 1D Gamma Array | %i x %i  |', size(gamma));

% Add GPU computation result
results{size(results,1)+1,1} = '5';
results{size(results,1),2} = 'Beam Profile 1D Gamma GPU Result within 1e-5';
results{size(results,1),3} = gpf;

% Add computation time
results{size(results,1)+1,1} = '6';
results{size(results,1),2} = 'Beam Profile 1D Gamma GPU Computation Time';
results{size(results,1),3} = gtime;

% Add CPU computation result
results{size(results,1)+1,1} = '7';
results{size(results,1),2} = 'Beam Profile 1D Gamma CPU Result within 1e-5';
results{size(results,1),3} = cpf;

% Add computation time
results{size(results,1)+1,1} = '8';
results{size(results,1),2} = 'Beam Profile 1D Gamma CPU Computation Time';
results{size(results,1),3} = ctime;

%% TEST 9/10/11/12: Simple 2D Gamma
%
% DESCRIPTION: This unit test creates a simple analytic reference dataset, 
%   modifies it by a known amount, and executes CalcGamma.  The gamma
%   result is then compared to a known solution
%
% RELEVANT REQUIREMENTS: none
%
% INPUT DATA: varargin{3}.gammaC
%
% CONDITION A (+): When computed using identical criteria, CalcGamma
%   produces an identical array to the reference with GPU Calculation
%
% CONDITION B (+): When computed using identical criteria, CalcGamma
%   produces an identical array to the reference with CPU Calculation
%
% CONDITION C (-): When computed using modified criteria, CalcGamma
%   produces a different result

% Declare Gamma criteria
percent = 3;
dta = 0.3;
local = 0;
res = 20;

% Declare size of data arrays
n = 100;

% Create reference data structure (peaks function)
ref.start = [0 0];
ref.width = [0.01 0.01];
ref.data = peaks(n) + 5;

% Declare modification values (to be applied to reference data)
shift = [0.05 0.05];
scale = 1.03;

% Modify reference data, saving to target structure
target.start = ref.start + shift;
target.width = ref.width;
target.data = ref.data * scale;

% Execute in try/catch statement
try

    % Execute CalcGamma using test dataset
    t = tic;
    gamma = CalcGamma(ref, target, percent, dta, 'local', local, ...
        'res', res);
    gtime = sprintf('%0.1f sec', toc(t));
    gpf = pass;
    cpf = pass;

    % If reference data exists
    if nargin == 3
        
        % If current value does not equal the reference
        if max(max(abs(gamma - varargin{3}.gammaC))) > 1e-5
            gpf = fail;
        end
        
        % Execute CalcGamma again using CPU
        t = tic;
        gamma = CalcGamma(ref, target, percent, dta, 'cpu', 1, ...
            'local', local, 'res', res);
        ctime = sprintf('%0.1f sec', toc(t));
        
        % If CPU value does not equal the reference
        if max(max(abs(gamma - varargin{3}.gammaC))) > 1e-5
            cpf = fail;
        end
        
        % Execute CalcGamma again, this time using modified criteria
        gamma = CalcGamma(ref, target, percent, dta * 2, ...
            'local', local, 'res', res);
        
        % If modified value equals the reference, the test fails
        if max(max(abs(gamma - varargin{3}.gammaC))) < 1e-5
            gpf = fail;
        end
        
        % Execute CalcGamma again, this time using modified criteria
        gamma = CalcGamma(ref, target, percent * 2, dta, ...
            'local', local, 'res', res);
        
        % If modified value equals the reference, the test fails
        if max(max(abs(gamma - varargin{3}.gammaC))) < 1e-5
            gpf = fail;
        end
        
    % Otherwise, set reference    
    elseif nargout == 4
        reference.gammaC = gamma;
        gpf = pass;
        cpf = pass;
    end
    
catch
    gpf = fail;
    cpf = fail;
end

% Add the dataset size to the preamble
preamble{length(preamble)+1} = ...
    sprintf('| Simple 2D Gamma Array | %i x %i  |', size(gamma));

% Add GPU computation result
results{size(results,1)+1,1} = '9';
results{size(results,1),2} = 'Simple 2D Gamma GPU Result within 1e-5';
results{size(results,1),3} = gpf;

% Add computation time
results{size(results,1)+1,1} = '10';
results{size(results,1),2} = 'Simple 2D Gamma GPU Computation Time';
results{size(results,1),3} = gtime;

% Add CPU computation result
results{size(results,1)+1,1} = '11';
results{size(results,1),2} = 'Simple 2D Gamma CPU Result within 1e-5';
results{size(results,1),3} = cpf;

% Add computation time
results{size(results,1)+1,1} = '12';
results{size(results,1),2} = 'Simple 2D Gamma CPU Computation Time';
results{size(results,1),3} = ctime;

%% TEST 13/14/15/16: Dose Volume 3D Gamma
%
% DESCRIPTION: This unit test creates a real world dose volume dataset 
%   and executes CalcGamma.  The gamma result is then compared to a known 
%   solution
%
% RELEVANT REQUIREMENTS: none
%
% INPUT DATA: varargin{2}.referenceD, varargin{2}.targetD, 
%   varargin{3}.gammaD
%
% CONDITION A (+): When computed using identical criteria, CalcGamma
%   produces an identical array to the reference with GPU Calculation
%
% CONDITION B (+): When computed using identical criteria, CalcGamma
%   produces an identical array to the reference with CPU Calculation

% Declare gamma criteria
percent = 3;
dta = 3;
local = 0;
restrict = 1;

% Load reference and target data from varargin
ref = varargin{2}.referenceD;
target = varargin{2}.targetD;

% Execute in try/catch statement
try

    % Execute CalcGamma using test dataset
    t = tic;
    gamma = CalcGamma(ref, target, percent, dta, 'local', local, ...
        'restrict', restrict);
    gtime = sprintf('%0.1f sec', toc(t));
    gpf = pass;
    cpf = pass;

    % If reference data exists
    if nargin == 3
        
        % If current value does not equal the reference
        if max(max(max(abs(gamma - varargin{3}.gammaD)))) > 1e-5
            gpf = fail;
        end
        
        % Execute CalcGamma again using CPU
        t = tic;
        gamma = CalcGamma(ref, target, percent, dta, 'cpu', 1, ...
            'local', local, 'restrict', restrict);
        ctime = sprintf('%0.1f sec', toc(t));
        
        % If CPU value does not equal the reference
        if max(max(max(abs(gamma - varargin{3}.gammaD)))) > 1e-5
            cpf = fail;
        end
        
    % Otherwise, set reference    
    elseif nargout == 4
        reference.gammaD = gamma;
        gpf = pass;
        cpf = pass;
    end
    
catch
    gpf = fail;
    cpf = fail;
end

% Add the dataset size to the preamble
preamble{length(preamble)+1} = ...
    sprintf('| Dose Volume 3D Gamma Array | %i x %i x %i |', size(gamma));

% Add GPU computation result
results{size(results,1)+1,1} = '13';
results{size(results,1),2} = 'Dose Volume 3D Gamma GPU Result within 1e-5';
results{size(results,1),3} = gpf;

% Add computation time
results{size(results,1)+1,1} = '14';
results{size(results,1),2} = 'Dose Volume 3D Gamma GPU Computation Time';
results{size(results,1),3} = gtime;

% Add CPU computation result
results{size(results,1)+1,1} = '15';
results{size(results,1),2} = 'Dose Volume 3D Gamma CPU Result within 1e-5';
results{size(results,1),3} = cpf;

% Add computation time
results{size(results,1)+1,1} = '16';
results{size(results,1),2} = 'Dose Volume 3D Gamma CPU Computation Time';
results{size(results,1),3} = ctime;

%% TEST 17: Data Validity
%
% DESCRIPTION: This test verifies that CalcGamma correctly fails when input
%   data is not valid.  This function tests error messages both with an
%   without the Event() function.
%
% RELEVANT REQUIREMENTS: none 
%
% INPUT DATA: No input data required
%
% CONDITION A (-): CalcGamma fails with too few arguments
%
% CONDITION B (-): CalcGamma fails if the reference structure is poorly
%   formatted
%
% CONDITION C (-): CalcGamma fails if the target structure is poorly
%   formatted
%
% CONDITION D (-): CalcGamma fails if the reference and target dimensions
%   differ
%
% CONDITION E (-): CalcGamma fails if the reference structure contains more
%   than three dimensions
%
% CONDITION F (-): CalcGamma fails if the target structure contains more
%   than three dimensions

% Start with passing result
pf = pass;

% Try to execute CalcGamma with too few arguments
try
    CalcGamma();
    pf = fail;
catch
    % The test passes if an error is thrown
end

% Initialize bad data
ref.data = rand(5);
ref.start = [0 0];
ref.width = [1 1 1];
target.data = rand(5);
target.start = [0 0];
target.width = [1 1 1];

% Try to execute CalcGamma with bad reference structure
try
    CalcGamma(ref, target, 1, 1);
    pf = fail;
catch
    % The test passes if an error is thrown
end

% Fix reference
ref.width = [1 1];

% Try to execute CalcGamma with bad target structure
try
    CalcGamma(ref, target, 1, 1);
    pf = fail;
catch
    % The test passes if an error is thrown
end

% Fix target
target.data = rand(5, 5, 5);
target.start = [0 0 0];
target.width = [1 1 1];

% Try to execute CalcGamma with different dimensions
try
    CalcGamma(ref, target, 1, 1);
    pf = fail;
catch
    % The test passes if an error is thrown
end

% Create 4-D datasets
ref.data = rand(5, 5, 5, 5);
ref.start = [0 0 0 0];
ref.width = [1 1 1 1];
target.data = rand(5, 5, 5, 5);
target.start = [0 0 0 0];
target.width = [1 1 1 1];

% Try to execute CalcGamma with different dimensions
try
    CalcGamma(ref, target, 1, 1);
    pf = fail;
catch
    % The test passes if an error is thrown
end

% Repeat without Event()
restoredefaultpath

% Try to execute CalcGamma with too few arguments
try
    CalcGamma();
    pf = fail;
catch
    % The test passes if an error is thrown
end

% Initialize bad data
ref.data = rand(5);
ref.start = [0 0];
ref.width = [1 1 1];
target.data = rand(5);
target.start = [0 0];
target.width = [1 1 1];

% Try to execute CalcGamma with bad reference structure
try
    CalcGamma(ref, target, 1, 1);
    pf = fail;
catch
    % The test passes if an error is thrown
end

% Fix reference
ref.width = [1 1];

% Try to execute CalcGamma with bad target structure
try
    CalcGamma(ref, target, 1, 1);
    pf = fail;
catch
    % The test passes if an error is thrown
end

% Fix target
target.data = rand(5, 5, 5);
target.start = [0 0 0];
target.width = [1 1 1];

% Try to execute CalcGamma with different dimensions
try
    CalcGamma(ref, target, 1, 1);
    pf = fail;
catch
    % The test passes if an error is thrown
end

% Create 4-D datasets
ref.data = rand(5, 5, 5, 5);
ref.start = [0 0 0 0];
ref.width = [1 1 1 1];
target.data = rand(5, 5, 5, 5);
target.start = [0 0 0 0];
target.width = [1 1 1 1];

% Try to execute CalcGamma with different dimensions
try
    CalcGamma(ref, target, 1, 1);
    pf = fail;
catch
    % The test passes if an error is thrown
end

% Add test result
results{size(results,1)+1,1} = '17';
results{size(results,1),2} = 'Data Validation Checks Functional';
results{size(results,1),3} = pf;

%% TEST 18/19: Code Analyzer Messages, Cyclomatic Complexity
%
% DESCRIPTION: This unit test uses the checkcode() MATLAB function to check
%   each function used by the application and return any Code Analyzer
%   messages that result.  The cumulative cyclomatic complexity is also
%   computed for each function and summed to determine the total
%   application complexity.  Although this test does not reference any
%   particular requirements, it is used during development to help identify
%   high risk code.
%
% RELEVANT REQUIREMENTS: none 
%
% INPUT DATA: No input data required
%
% CONDITION A (+): Report any code analyzer messages for all functions
%   called by CalcGamma
%
% CONDITION B (+): Report the cumulative cyclomatic complexity for all
%   functions called by CalcGamma

% Search for required functions
fList = matlab.codetools.requiredFilesAndProducts('CalcGamma.m');

% Initialize complexity and messages counters
comp = 0;
mess = 0;

% Loop through each dependency
for i = 1:length(fList)
    
    % Execute checkcode
    inform = checkcode(fList{i}, '-cyc');
    
    % Loop through results
    for j = 1:length(inform)
       
        % Check for McCabe complexity output
        c = regexp(inform(j).message, ...
            '^The McCabe complexity .+ is ([0-9]+)\.$', 'tokens');
        
        % If regular expression was found
        if ~isempty(c)
            
            % Add complexity
            comp = comp + str2double(c{1});
            
        else
            
            % If not an invalid code message
            if ~strncmp(inform(j).message, 'Filename', 8)
                
                % Log message
                Event(sprintf('%s in %s', inform(j).message, fList{i}), ...
                    'CHCK');

                % Add as code analyzer message
                mess = mess + 1;
            end
        end
        
    end
end

% Add code analyzer messages counter to results
results{size(results,1)+1,1} = '18';
results{size(results,1),2} = 'Code Analyzer Messages';
results{size(results,1),3} = sprintf('%i', mess);

% Add complexity results
results{size(results,1)+1,1} = '19';
results{size(results,1),2} = 'Cumulative Cyclomatic Complexity';
results{size(results,1),3} = sprintf('%i', comp);

%% Finish up
% Close all figures
close all force;

% If no return variables are present, print the results
if nargout == 0
    
    % Print preamble
    for j = 1:length(preamble)
        fprintf('%s\n', preamble{j});
    end
    fprintf('\n');
    
    % Loop through each table row
    for j = 1:size(results,1)
        
        % Print table row
        fprintf('| %s |\n', strjoin(results(j,:), ' | '));
       
        % If this is the first column
        if j == 1
            
            % Also print a separator row
            fprintf('|%s\n', repmat('----|', 1, size(results,2)));
        end

    end
    fprintf('\n');
    
    % Print footnotes
    for j = 1:length(footnotes) 
        fprintf('%s<br>\n', footnotes{j});
    end
    
% Otherwise, return the results as variables    
else

    % Store return variables
    if nargout >= 1; varargout{1} = preamble; end
    if nargout >= 2; varargout{2} = results; end
    if nargout >= 3; varargout{3} = footnotes; end
    if nargout >= 4; varargout{4} = reference; end
end
