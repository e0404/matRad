function outDose = matRad_readOpenStack(folder)

if contains(folder,'*')% && all(modulation ~= false)
    folder = dir(folder);
    for i = 1:length(folder)
        folders{i} = [folder(i).folder filesep folder(i).name];
    end
else
    folders{1} = folder;
end


for f = 1:length(folders)
    currFolder = folders{f};
    topasConfig = MatRad_TopasConfig();
    topasCubes = matRad_readTopasData(currFolder);
    
    load([currFolder filesep 'dij.mat']);
    load([currFolder filesep 'weights.mat']);
    
    ctScen = 1;
    fnames = fieldnames(topasCubes);
    dij.MC_tallies = fnames;
    
    if ~isfield(topasCubes,'RBE')
        for f = 1:numel(fnames)
            dij.(fnames{f}){ctScen,1} = sum(w(:,ctScen))*reshape(topasCubes.(fnames{f}),[],1);
        end
    else
        for d = 1:length(stf)
            dij.physicalDose{ctScen,1}(:,d)    = sum(w)*reshape(topasCubes.(['physicalDose_beam',num2str(d)]),[],1);
            dij.alpha{ctScen,1}(:,d)           = reshape(topasCubes.(['alpha_beam',num2str(d)]),[],1);
            dij.beta{ctScen,1}(:,d)            = reshape(topasCubes.(['beta_beam',num2str(d)]),[],1);
            
            [dij.ax,dij.bx] = matRad_getPhotonLQMParameters(cst,prod(ct.cubeDim),1);
            dij.abx(dij.bx>0) = dij.ax(dij.bx>0)./dij.bx(dij.bx>0);
            
            dij.mAlphaDose{ctScen,1}(:,d)      = dij.physicalDose{ctScen,1}(:,d) .* dij.alpha{ctScen,1}(:,d);
            dij.mSqrtBetaDose{ctScen,1}(:,d)   = sqrt(dij.physicalDose{ctScen,1}(:,d)) .* dij.beta{ctScen,1}(:,d);
        end
    end
    
    resultGUI    = matRad_calcCubes(w,dij,1);
    
    if ~exist('outDose','var')
        outDose = zeros(size(resultGUI.physicalDose));
    end
        
    outDose = outDose + resultGUI.physicalDose/length(folders);
end