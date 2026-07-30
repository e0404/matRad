function test_suite = test_roundCompat
% The output should always be test_suite, and the function name the same as
% your file name

% To collect all tests defined below, this is needed in newer Matlab
% versions. test_functions will collect function handles to below test
% functions
test_functions = localfunctions();

% This will initialize the test suite, i.e., take the functions from
% test_functions, check if they contain "test", convert them into a MOxUnit
% Test Case, and add them to the test-runner
initTestSuite;

function test_roundCompatWithoutDigits
% no digit count given behaves like plain round
assertEqual(matRad_roundCompat(2.4), 2);
assertEqual(matRad_roundCompat(2.6), 3);
assertEqual(matRad_roundCompat(-2.6), -3);
assertEqual(matRad_roundCompat([1.2, 1.8; -0.4, 2.5]), [1, 2; 0, 3]);
% an empty digit count is the same as omitting it
assertEqual(matRad_roundCompat(2.6, []), 3);

function test_roundCompatPositiveDigits
assertElementsAlmostEqual(matRad_roundCompat(3.14159, 2), 3.14, 'absolute', 1e-12);
assertElementsAlmostEqual(matRad_roundCompat(3.14159, 4), 3.1416, 'absolute', 1e-12);
% zero digits rounds to integers
assertEqual(matRad_roundCompat(2.6, 0), 3);

function test_roundCompatNegativeDigits
% negative digit counts round to the left of the decimal point, which
% matRad_StfGeneratorParticleIMPT-style range shifter placement relies on
assertElementsAlmostEqual(matRad_roundCompat(123.4, -1), 120, 'absolute', 1e-12);
assertElementsAlmostEqual(matRad_roundCompat(1234, -2), 1200, 'absolute', 1e-12);

function test_roundCompatVectorized
angles = [-180, -172.5, 0, 12.34567, 179.99999];
expected = [-180, -172.5, 0, 12.3457, 180];
assertElementsAlmostEqual(matRad_roundCompat(angles, 4), expected, 'absolute', 1e-9);
% shape is preserved
assertEqual(size(matRad_roundCompat(reshape(1:6, 2, 3) + 0.123456, 3)), [2, 3]);

function test_roundCompatMatchesNativeRoundOnMatlab
% where Matlab's two-argument round exists, the compat function must agree
% with it, so that results do not depend on which one a call site uses
matRad_cfg = MatRad_Config.instance();
if matRad_cfg.isOctave
    moxunit_throw_test_skipped_exception('Two-argument round does not exist in Octave');
end

values = [pi, -pi, 0, 1e-6, 12345.6789, 0.5, -0.5, 2.675];
for n = -2:6
    assertEqual(matRad_roundCompat(values, n), round(values, n), ...
                sprintf('mismatch for n = %d', n));
end

function test_roundCompatSingleArgumentRoundIsOctaveSafe
% guard against reintroducing round(x,n): every call outside the compat
% function itself must go through matRad_roundCompat, since Octave's round
% rejects the digit argument
matRad_cfg = MatRad_Config.instance();
offenders = {};
% tests and examples are run under Octave as well, and the call that
% prompted this function lived in a test
roots = {matRad_cfg.matRadSrcRoot, ...
         fullfile(matRad_cfg.matRadRoot, 'test'), ...
         matRad_cfg.exampleFolder};
files = {};
for r = 1:numel(roots)
    files = [files, helper_gatherMFiles(roots{r})]; %#ok<AGROW>
end

for i = 1:numel(files)
    % this file necessarily names the pattern it forbids
    if ~isempty(strfind(files{i}, 'roundCompat'))
        continue
    end
    lines = helper_readLines(files{i});
    for j = 1:numel(lines)
        if helper_hasTwoArgumentRound(lines{j})
            offenders{end + 1} = sprintf('%s:%d', files{i}, j); %#ok<AGROW>
        end
    end
end
assertTrue(isempty(offenders), ...
           sprintf('round(x,n) is not available in Octave, use matRad_roundCompat:\n%s', ...
                   strjoin(offenders, sprintf('\n'))));

function files = helper_gatherMFiles(folder)
files = {};
entries = dir(folder);
for i = 1:numel(entries)
    name = entries(i).name;
    if strcmp(name, '.') || strcmp(name, '..')
        continue
    end
    full = fullfile(folder, name);
    if entries(i).isdir
        files = [files, helper_gatherMFiles(full)]; %#ok<AGROW>
    elseif numel(name) > 2 && strcmp(name(end - 1:end), '.m')
        files{end + 1} = full; %#ok<AGROW>
    end
end

function lines = helper_readLines(file)
fid = fopen(file, 'r');
if fid < 0
    lines = {};
    return
end
raw = fread(fid, inf, '*char')';
fclose(fid);
lines = strsplit(raw, sprintf('\n'));

function tf = helper_hasTwoArgumentRound(line)
% true if the line calls round() with more than one top level argument,
% ignoring comments
tf = false;
trimmed = strtrim(line);
if isempty(trimmed) || strcmp(trimmed(1), '%')
    return
end

starts = regexp(line, '(?<![\w.])round\s*\(', 'end');
for s = 1:numel(starts)
    depth = 1;
    k = starts(s) + 1;
    while k <= numel(line) && depth > 0
        switch line(k)
            case {'(', '[', '{'}
                depth = depth + 1;
            case {')', ']', '}'}
                depth = depth - 1;
            case ','
                if depth == 1
                    tf = true;
                    return
                end
        end
        k = k + 1;
    end
end
