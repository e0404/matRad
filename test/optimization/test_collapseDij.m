function test_suite = test_collapseDij
    
    test_functions=localfunctions();
    
    initTestSuite;
    
function test_collapse_dij_numeric_photon
    % matrix
    testData = load('photons_testData.mat');

    dijNew = matRad_collapseDij(testData.dij);
    assertEqual(dijNew.totalNumOfBixels, testData.dij.numOfBeams);
    assertEqual(dijNew.totalNumOfRays, testData.dij.numOfBeams);
    assertEqual(dijNew.numOfBeams, testData.dij.numOfBeams);
    assertEqual(dijNew.numOfRaysPerBeam, ones(testData.dij.numOfBeams,1));
    assertEqual(dijNew.beamNum, transpose(1:testData.dij.numOfBeams));
    assertEqual(dijNew.bixelNum, ones(testData.dij.numOfBeams, 1));
    assertEqual(dijNew.rayNum, ones(testData.dij.numOfBeams, 1));
    assertEqual(dijNew.doseGrid, testData.dij.doseGrid);
    assertEqual(dijNew.ctGrid, testData.dij.ctGrid);

    assertEqual(size(dijNew.physicalDose), size(testData.dij.physicalDose));
    assertEqual(cellfun(@isempty, dijNew.physicalDose), cellfun(@isempty, testData.dij.physicalDose));
    
    for i = 1:numel(dijNew.physicalDose)
        if ~isempty(dijNew.physicalDose{i})
            assertEqual(size(dijNew.physicalDose{i},1), size(testData.dij.physicalDose{i},1));
            assertEqual(size(dijNew.physicalDose{i},2), testData.dij.numOfBeams);
            assertElementsAlmostEqual(sum(dijNew.physicalDose{i}(:)), sum(testData.dij.physicalDose{i}(:)), 'relative', 1e-5, 1e-5);
        end
    end

function test_collapse_dij_numeric_proton
    % matrix
    testData = load('protons_testData.mat');

    dijNew = matRad_collapseDij(testData.dij);
    assertEqual(dijNew.totalNumOfBixels, testData.dij.numOfBeams);
    assertEqual(dijNew.totalNumOfRays, testData.dij.numOfBeams);
    assertEqual(dijNew.numOfBeams, testData.dij.numOfBeams);
    assertEqual(dijNew.numOfRaysPerBeam, ones(testData.dij.numOfBeams,1));
    assertEqual(dijNew.beamNum, transpose(1:testData.dij.numOfBeams));
    assertEqual(dijNew.bixelNum, ones(testData.dij.numOfBeams, 1));
    assertEqual(dijNew.rayNum, ones(testData.dij.numOfBeams, 1));
    assertEqual(numel(dijNew.numParticlesPerMU),testData.dij.numOfBeams);
    assertEqual(sum(dijNew.numParticlesPerMU),sum(testData.dij.numParticlesPerMU));
    assertEqual(dijNew.doseGrid, testData.dij.doseGrid);
    assertEqual(dijNew.ctGrid, testData.dij.ctGrid);

    assertEqual(size(dijNew.physicalDose), size(testData.dij.physicalDose));
    assertEqual(cellfun(@isempty, dijNew.physicalDose), cellfun(@isempty, testData.dij.physicalDose));    
    
    quantities = {'physicalDose','mLETDose'};

    for q = 1:numel(quantities)    
        for i = 1:numel(dijNew.(quantities{q}))
            if ~isempty(dijNew.(quantities{q}))
                assertEqual(size(dijNew.(quantities{q}){i},1), size(testData.dij.(quantities{q}){i},1));
                assertEqual(size(dijNew.(quantities{q}){i},2), testData.dij.numOfBeams);
                assertElementsAlmostEqual(sum(dijNew.(quantities{q}){i}(:)), sum(testData.dij.(quantities{q}){i}(:)), 'relative', 1e-5, 1e-5);
            end
        end
    end

function test_collapse_dij_numeric_carbon
    % matrix
    testData = load('carbon_testData.mat');

    dijNew = matRad_collapseDij(testData.dij);
    assertEqual(dijNew.totalNumOfBixels, testData.dij.numOfBeams);
    assertEqual(dijNew.totalNumOfRays, testData.dij.numOfBeams);
    assertEqual(dijNew.numOfBeams, testData.dij.numOfBeams);
    assertEqual(dijNew.numOfRaysPerBeam, ones(testData.dij.numOfBeams,1));
    assertEqual(dijNew.beamNum, transpose(1:testData.dij.numOfBeams));
    assertEqual(dijNew.bixelNum, ones(testData.dij.numOfBeams, 1));
    assertEqual(dijNew.rayNum, ones(testData.dij.numOfBeams, 1));
    assertEqual(numel(dijNew.numParticlesPerMU),testData.dij.numOfBeams);
    assertEqual(sum(dijNew.numParticlesPerMU),sum(testData.dij.numParticlesPerMU));
    assertEqual(dijNew.doseGrid, testData.dij.doseGrid);
    assertEqual(dijNew.ctGrid, testData.dij.ctGrid);

    assertEqual(size(dijNew.physicalDose), size(testData.dij.physicalDose));
    assertEqual(cellfun(@isempty, dijNew.physicalDose), cellfun(@isempty, testData.dij.physicalDose));    
    
    quantities = {'physicalDose','mLETDose','mAlphaDose','mSqrtBetaDose'};

    for q = 1:numel(quantities)    
        for i = 1:numel(dijNew.(quantities{q}))
            if ~isempty(dijNew.(quantities{q}))
                assertEqual(size(dijNew.(quantities{q}){i},1), size(testData.dij.(quantities{q}){i},1));
                assertEqual(size(dijNew.(quantities{q}){i},2), testData.dij.numOfBeams);
                assertElementsAlmostEqual(sum(dijNew.(quantities{q}){i}(:)), sum(testData.dij.(quantities{q}){i}(:)), 'relative', 1e-5, 1e-5);
            end
        end
    end