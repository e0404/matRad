function resultGUI = matRad_readOpenStack(folder)

if contains(folder,'*')% && all(modulation ~= false)
    folder = dir(folder);
    for i = 1:length(folder)
        if any(isstrprop(folder(i).name,'digit'))
            folders{i} = [folder(i).folder filesep folder(i).name];
        end
    end
else
    folders{1} = folder;
end
folders = folders(~cellfun('isempty',folders));

for f = 1:length(folders)
    try
        currFolder = folders{f};
        topasConfig = MatRad_TopasConfig();
        topasCubes = matRad_readTopasData(currFolder);
        
        loadedVars = load([currFolder filesep 'dij.mat'],'dij');
        dij = loadedVars.dij;
        loadedVars = load([currFolder filesep 'weights.mat'],'w');
        w = loadedVars.w;
        
        ctScen = 1;
        fnames = fieldnames(topasCubes);
        dij.MC_tallies = fnames;
        
        if ~isfield(topasCubes,'alpha_beam1')
            for f = 1:numel(fnames)
                dij.(fnames{f}){ctScen,1} = sum(w(:,ctScen))*reshape(topasCubes.(fnames{f}),[],1);
            end
        else
            for d = 1:dij.numOfBeams
                dij.physicalDose{ctScen,1}(:,d)    = sum(w)*reshape(topasCubes.(['physicalDose_beam',num2str(d)]),[],1);
                dij.alpha{ctScen,1}(:,d)           = reshape(topasCubes.(['alpha_beam',num2str(d)]),[],1);
                dij.beta{ctScen,1}(:,d)            = reshape(topasCubes.(['beta_beam',num2str(d)]),[],1);
                
                dij.mAlphaDose{ctScen,1}(:,d)      = dij.physicalDose{ctScen,1}(:,d) .* dij.alpha{ctScen,1}(:,d);
                dij.mSqrtBetaDose{ctScen,1}(:,d)   = sqrt(dij.physicalDose{ctScen,1}(:,d)) .* dij.beta{ctScen,1}(:,d);
            end
        end
        
        if length(folders) > 1
            outDose    = matRad_calcCubes(ones(dij.numOfBeams,1),dij,1);
            if ~exist('resultGUI')
                for i = 1:dij.numOfBeams
                    beamInfo(i).suffix = ['_beam', num2str(i)];
                end
                beamInfo(dij.numOfBeams+1).suffix = '';
                for i = 1:length(beamInfo)
                    resultGUI.(['physicalDose', beamInfo(i).suffix]) = zeros(dij.ctGrid.dimensions);
                    resultGUI.(['RBExD', beamInfo(i).suffix]) = zeros(dij.ctGrid.dimensions);
                    
                    resultGUI.(['alpha', beamInfo(i).suffix]) = {};
                    resultGUI.(['beta', beamInfo(i).suffix]) = {};
                    resultGUI.(['RBE', beamInfo(i).suffix]) = {};
                    resultGUI.(['effect', beamInfo(i).suffix]) = {};
                end
            end
            
            for i = 1:length(beamInfo)
                resultGUI.(['physicalDose', beamInfo(i).suffix]) = resultGUI.(['physicalDose', beamInfo(i).suffix]) + outDose.(['physicalDose', beamInfo(i).suffix])/length(folders);
                resultGUI.(['RBExD', beamInfo(i).suffix]) = resultGUI.(['RBExD', beamInfo(i).suffix]) + outDose.(['RBExD', beamInfo(i).suffix])/length(folders);
                
                resultGUI.(['alpha', beamInfo(i).suffix]){f} = outDose.(['alpha', beamInfo(i).suffix]);
                resultGUI.(['beta', beamInfo(i).suffix]){f} = outDose.(['beta', beamInfo(i).suffix]);
                resultGUI.(['RBE', beamInfo(i).suffix]){f} = outDose.(['RBE', beamInfo(i).suffix]);
                resultGUI.(['effect', beamInfo(i).suffix]){f} = outDose.(['effect', beamInfo(i).suffix]);
                resultGUI.samples = f;
            end
        else
            resultGUI    = matRad_calcCubes(ones(dij.numOfBeams,1),dij,1);
        end
    catch
        warning(['error in file ',currFolder]);
    end
end

end