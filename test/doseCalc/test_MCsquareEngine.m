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
        % linesConfigFile = readlines(fullfile(matRad_cfg.primaryUserFolder, 'MCsquare', 'MCsquareConfig.txt'));
        fid = fopen(fullfile(matRad_cfg.primaryUserFolder, 'MCsquare', 'MCsquareConfig.txt'),'r');
            linesConfigFile = {};
            while ~feof(fid)
                linesConfigFile{end+1,1} = fgetl(fid);
            end
        fclose(fid);


        assertTrue(any(contains(linesConfigFile, "Num_Primaries 42")));
        
        % Read currBixel file
        % linesBixelFile = readlines(fullfile(matRad_cfg.primaryUserFolder, 'MCsquare', 'currBixels.txt'));
        fid = fopen(fullfile(matRad_cfg.primaryUserFolder, 'MCsquare', 'currBixels.txt'),'r');
            linesBixelFile = {};
            while ~feof(fid)
                linesBixelFile{end+1,1} = fgetl(fid);
            end
        fclose(fid);

        assertTrue(any(contains(linesBixelFile, "##NumberOfFields")));
        assertTrue(str2double(linesBixelFile(find(strcmp(linesBixelFile, "##NumberOfFields"))+1)) == numel(stf));

    end


