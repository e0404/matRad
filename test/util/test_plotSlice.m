function test_suite = test_plotSlice

test_functions=localfunctions();

initTestSuite;

function test_plot_ct_only

    load BOXPHANTOM.mat
    [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSlice(ct);
    assertTrue(isempty(hCMap));
    assertTrue(isempty(hDose));
    assertFalse(isempty(hCt));
    if ~moxunit_util_platform_is_octave
       assertTrue(isa(hCt, 'matlab.graphics.primitive.Image'))
    end
    assertTrue(isempty(hContour));
    assertTrue(isempty(hIsoDose));
    close(gcf);

    load PROSTATE.mat
    [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSlice(ct, 'slice', 91, 'cst', cst);
    assertTrue(isempty(hCMap));
    assertTrue(isempty(hDose));
    assertFalse(isempty(hCt));
    if ~moxunit_util_platform_is_octave
       assertTrue(isa(hCt, 'matlab.graphics.primitive.Image'))
    end
    assertTrue(isa(hContour, 'cell'));
    assertTrue(isempty(hIsoDose));
    close(gcf);

function test_plot_dose_slice

    load protons_testData.mat
    [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSlice(ct, 'dose', resultGUI.physicalDose);
    assertFalse(isempty(hCt));
    assertTrue(isempty(hContour));
    assertTrue(isempty(hIsoDose));
    if ~moxunit_util_platform_is_octave
        assertTrue(isa(hCMap, 'matlab.graphics.illustration.ColorBar'));
        assertTrue(isa(hDose, 'matlab.graphics.primitive.Image'));
        assertTrue(isa(hCt, 'matlab.graphics.primitive.Image'))
    end
    close(gcf);

    load helium_testData.mat
    [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSlice(ct, 'dose', resultGUI.physicalDose);
    assertFalse(isempty(hCt));
    assertTrue(isempty(hContour));
    assertTrue(isempty(hIsoDose));
    if ~moxunit_util_platform_is_octave
        assertTrue(isa(hCMap, 'matlab.graphics.illustration.ColorBar'));
        assertTrue(isa(hDose, 'matlab.graphics.primitive.Image'));
        assertTrue(isa(hCt, 'matlab.graphics.primitive.Image'))
    end
    close(gcf);

    load carbon_testData.mat
    [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSlice(ct, 'dose', resultGUI.physicalDose);
    assertFalse(isempty(hCt));
    assertTrue(isempty(hContour));
    assertTrue(isempty(hIsoDose));
    if ~moxunit_util_platform_is_octave
        assertTrue(isa(hCMap, 'matlab.graphics.illustration.ColorBar'));
        assertTrue(isa(hDose, 'matlab.graphics.primitive.Image'));
        assertTrue(isa(hCt, 'matlab.graphics.primitive.Image'))
    end
    close(gcf);

    load photons_testData.mat
    [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSlice(ct, 'dose', resultGUI.physicalDose, 'plane', 3);
    assertFalse(isempty(hCt));
    assertTrue(isempty(hContour));
    assertTrue(isempty(hIsoDose));
    if ~moxunit_util_platform_is_octave
        assertTrue(isa(hCMap, 'matlab.graphics.illustration.ColorBar'));
        assertTrue(isa(hDose, 'matlab.graphics.primitive.Image'));
        assertTrue(isa(hCt, 'matlab.graphics.primitive.Image'))
    end
    close(gcf);

function test_optional_input
    
    load photons_testData.mat
    hF = figure();
    doseCube = resultGUI.physicalDose;
    boolVOIselection = ones(1, size(cst, 1));
    [hCMap,hDose,hCt,hContour,hIsoDose] = matRad_plotSlice(ct,  ...
        'dose', doseCube, 'axesHandle', gca,                    ...
        'cst', cst, 'cubeIdx', 1, 'plane', 3, 'slice', 5,      ...
        'thresh', 0.1*max(doseCube(:)), 'alpha', 0.8,           ...
        'contourColorMap', white, 'doseColorMap', jet,          ...
        'doseWindow', [min(doseCube(:)) 1.1*max(doseCube(:))],  ...
        'doseIsoLevels', max(doseCube(:)).*[0.5, 0.6, 0.7, 0.8, 0.9, 0.92, 0.95, 0.97, 0.99], ...
        'voiSelection', boolVOIselection,                       ...
        'colorBarLabel', 'Absorbed Dose [Gy]',                  ...
        'boolPlotLegend', 1, 'showCt', 0, 'FontSize', 13);

    assertTrue(isempty(hCt));
    assertTrue(isa(hContour, "cell"));
    assertTrue(isa(hIsoDose, "cell"));
    if ~moxunit_util_platform_is_octave
        assertTrue(isa(hDose, 'matlab.graphics.primitive.Image'));
    end
    close(hF);

function test_title_input
    load photons_testData.mat
    hF = figure();
    hA = axes(hF);
    tString = 'Hello';
    matRad_plotSlice(ct,'title','Hello','axesHandle',hA);
    assertTrue(isequal('Hello',get(get(hA,'Title'),'String')));
    close(hF);
