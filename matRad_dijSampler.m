function [cst,dij] = matRad_dijSampler(cst,ct,dij, margin_size,sparCity)
% marginsize given in pixel thickness
% fraction of voxels to be sampled [ Target OAR]

m = (2*margin_size)+1;
msk = ones(m,m,m);    % mask for convolution

volThresh = 1000;  % setting full sampling for small structures 
idx = [];
fullmap =[];            % just for displaying 
s = RandStream('mlfg6331_64');                               % initialize random number generator (seed value )
cst = matRad_setOverlapPriorities(cst,ct.cubeDim);

for i= 1:size(cst,1) 
    if  ~isempty(cst{i,6})
        for j=1:numel(cst{i,6})
            cst{i,6}(j).numOfVoxels = numel(cst{i,4}{1});
        end 
    end
end

for i=1:size(cst,1)
    
    % fraction of voxels sampled 
    if strcmp(cst{i,3}, 'TARGET') || (numel(cst{i,4}{1}) < volThresh)                % threshold volume for small structures
        sp = sparCity(1);
    elseif strcmp(cst{i,3}, 'OAR')
        sp = sparCity(2);                                   % sampling fraction 
    end
    % getting the margins 
    mVOI = zeros(ct.cubeDim) ;
    mVOI( cst{i,4}{1}) = 1;
    mVc = convn(mVOI,msk,'same');
    mVd = mVc - mVOI.*sum(msk(:));   % finding closed edge voxels
    idx = find(mVd<0);   % idx
    % randomly sampling within the margins
    mVOI(idx) = 0;
    list = find(mVOI);
    idx2 = datasample(s, list, floor(sp * numel(list)),'Replace',false );
    cst{i,4}{2} = cst{i,4}{1};                                                      % saving the original idx ... required for ray casting !!!
    cst{i,4}{1} = [idx; idx2];
   
    idxRm =  setdiff(cst{i,4}{2},cst{i,4}{1});
    
    fullmap = [ fullmap; idx; idx2];                                    % for plotting the map of selected voxels
    
end
 vec = zeros( numel(mVOI),1);
 vec (fullmap)  = 1;
tic,  dij.physicalDose{1} = ( dij.physicalDose{1} .* vec);toc

     if isfield(dij,'mAlphaDose')        
        dij.mAlphaDose{1} = dij.mAlphaDose{1}.* vec;
        dij.mSqrtBetaDose{1} = dij.mSqrtBetaDose{1}.* vec;        
    end
xx = zeros(ct.cubeDim);
xx(fullmap) =1;
figure, imagesc(xx(:,:,80));
end


