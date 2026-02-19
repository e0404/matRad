function test_suite = test_MCsquareEngine

test_functions=localfunctions();

initTestSuite;


function test_MCsquareDoseCalcBasic

    matRad_cfg = MatRad_Config.instance();
    radModes = DoseEngines.matRad_ParticleMCsquareEngine.possibleRadiationModes;

    for i = 1:numel(radModes)
        load([radModes{i} '_testData.mat']);
        pln.bioModel = matRad_bioModel(radModes{i},'none');
        
        w = ones(1,sum([stf(:).totalNumOfBixels]));
    
        pln.propDoseCalc.engine = 'MCsquare';
        pln.propDoseCalc.externalCalculation = 'write';
        pln.propDoseCalc.numHistoriesDirect = 42;
        resultGUI = matRad_calcDoseForward(ct,cst,stf,pln, w);

        assertTrue(exist(fullfile(matRad_cfg.primaryUserFolder, 'MCsquare'), 'dir')==7); % Check it exists and its a folder
        assertTrue(exist(fullfile(matRad_cfg.primaryUserFolder, 'MCsquare', 'MCsquareConfig.txt'), 'file')==2);
        assertTrue(exist(fullfile(matRad_cfg.primaryUserFolder, 'MCsquare', 'currBixels.txt'), 'file')==2);
        
        
        % Check parameters
        % Read config file
        fid = fopen(fullfile(matRad_cfg.primaryUserFolder, 'MCsquare', 'MCsquareConfig.txt'),'r');
            linesConfigFile = {};
            while ~feof(fid)
                linesConfigFile{end+1,1} = fgetl(fid);
            end
        fclose(fid);


        assertTrue(any(strcmp(linesConfigFile, "Num_Primaries 42")));
        
        % Read currBixel file
        fid = fopen(fullfile(matRad_cfg.primaryUserFolder, 'MCsquare', 'currBixels.txt'),'r');
            linesBixelFile = {};
            while ~feof(fid)
                linesBixelFile{end+1,1} = fgetl(fid);
            end
        fclose(fid);

        assertTrue(any(strcmp(linesBixelFile, "##NumberOfFields")));
        assertTrue(str2double(linesBixelFile(find(strcmp(linesBixelFile, "##NumberOfFields"))+1)) == numel(stf));

    end

function test_MCsquareRaShi
    matRad_cfg = MatRad_Config.instance();
    radModes = DoseEngines.matRad_ParticleMCsquareEngine.possibleRadiationModes;

    for i = 1:numel(radModes)
        load([radModes{i} '_testData.mat']);
        pln.bioModel = matRad_bioModel(radModes{i},'none');

        pln.propStf.gantryAngles = 0;
        pln.propStf.couchAngles = 0;
        stfGenerator = matRad_StfGeneratorParticleSingleBeamlet(pln);
        stfGenerator.useRangeShifter = true;
        stf = stfGenerator.generate(ct,cst);


        w = 1;

        pln.propDoseCalc.engine = 'MCsquare';
        pln.propDoseCalc.externalCalculation = 'write';
        pln.propDoseCalc.numHistoriesDirect = 42;
        resultGUI = matRad_calcDoseForward(ct,cst,stf,pln, w);

        assertTrue(exist(fullfile(matRad_cfg.primaryUserFolder, 'MCsquare'), 'dir')==7); % Check it exists and its a folder
        assertTrue(exist(fullfile(matRad_cfg.primaryUserFolder, 'MCsquare', 'MCsquareConfig.txt'), 'file')==2);
        assertTrue(exist(fullfile(matRad_cfg.primaryUserFolder, 'MCsquare', 'currBixels.txt'), 'file')==2);
        
        fid = fopen(fullfile(matRad_cfg.primaryUserFolder, 'MCsquare', 'currBixels.txt'),'r');
            linesBixelFile = {};
            while ~feof(fid)
                linesBixelFile{end+1,1} = fgetl(fid);
            end
        fclose(fid);

        assertTrue(any(strcmp(linesBixelFile, "####RangeShifterWaterEquivalentThickness")));
        assertTrue(str2double(linesBixelFile(find(strcmp(linesBixelFile, "####RangeShifterWaterEquivalentThickness"))+1)) == stf.ray.rangeShifter.eqThickness);
    end