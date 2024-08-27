function test_suite = test_themes
%The output should always be test_suite, and the function name the same as
%your file name

%% Header
% The header is required to be in this format for automatic test collection
% by MOxUnit

%To collect all tests defined below, this is needed in newer Matlab
%versions. test_functions will collect function handles to below test
%functions
test_functions=localfunctions(); 

% This will initialize the test suite, i.e., take the functions from
% test_functions, check if they contain "test", convert them into a MOxUnit
% Test Case, and add them to the test-runner
initTestSuite;

%% Custom Tests
% Tests use assert*-like Functions to check outputs etc:
% assertTrue(a) - a is true
% assertFalse(a) - a is false
% assertEqual(a,b) - a and be are equal (isequal)
% assertElementsAlmostEqual(a,b) - numerical test for all vector / matrix elements. Has Additional arguments for absolute / relative tolerance 
% assertVectorsAlmostEqual(a,b) - numerical test using vector norm
% assertExceptionThrown(f,id) - test if exception of id is thrown. Take care of Octave issues with exception id (or don't provide id)
% Check MOxUnit for more information or look at other tests

function themes = helper_getThemes()
    themes = {?matRad_ThemeDark, ?matRad_ThemeLight};

% Test constructor
function test_light_constructor
    themes = helper_getThemes();
    for i = 1:length(themes)
        theme = str2func(themes{i}.Name);
        theme = theme();
        assertTrue(isa(theme, 'matRad_Theme'));
    end

function test_properties_struct
    themes = helper_getThemes();
    for i = 1:length(themes)
        themeC = str2func(themes{i}.Name);
        themeC = themeC();
        theme = struct(themeC);            

        % Verify theme name is not empty
        assertFalse(isempty(theme.name));
        assertTrue(ischar(theme.name));
        assertTrue(isrow(theme.name));
        assertEqual(theme.name, themeC.name);

        % Verify background color size and range
        assertEqual(size(theme.backgroundColor), [1 3]);
        assertTrue(all(theme.backgroundColor >= 0));
        assertTrue(all(theme.backgroundColor <= 1));
        assertEqual(theme.backgroundColor, themeC.backgroundColor);

        % Verify element color size and range
        assertEqual(size(theme.elementColor), [1 3]);
        assertTrue(all(theme.elementColor >= 0));
        assertTrue(all(theme.elementColor <= 1));
        assertEqual(theme.elementColor, themeC.elementColor);

        % Verify text color size and range
        assertEqual(size(theme.textColor), [1 3]);
        assertTrue(all(theme.textColor >= 0));
        assertTrue(all(theme.textColor <= 1));
        assertEqual(theme.textColor, themeC.textColor);

        % Verify highlight color size and range
        assertEqual(size(theme.highlightColor), [1 3]);
        assertTrue(all(theme.highlightColor >= 0));
        assertTrue(all(theme.highlightColor <= 1));
        assertEqual(theme.highlightColor, themeC.highlightColor);

        % Verify font size is a positive scalar
        assertTrue(isscalar(theme.fontSize));
        assertTrue(isnumeric(theme.fontSize));
        assertTrue(theme.fontSize > 0);
        assertEqual(theme.fontSize, themeC.fontSize);

        % Verify font weight is not empty and either 'normal' or 'bold'
        assertFalse(isempty(theme.fontWeight));            
        assertTrue(any(strcmp(theme.fontWeight, {'normal', 'bold'})));
        assertEqual(theme.fontWeight, themeC.fontWeight);

        % Verify font name is not empty
        assertFalse(isempty(theme.fontName));
        assertTrue(ischar(theme.fontName));
        assertTrue(isrow(theme.fontName));
        assertEqual(theme.fontName, themeC.fontName);

        % Verify author is not empty
        assertFalse(isempty(theme.author));
        assertTrue(ischar(theme.author));
        assertTrue(isrow(theme.author));
        assertEqual(theme.author, themeC.author);

        % Verify description is not empty
        assertFalse(isempty(theme.description));
        assertTrue(ischar(theme.description));
        assertTrue(isrow(theme.description));
        assertEqual(theme.description, themeC.description);
    end
