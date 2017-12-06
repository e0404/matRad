function matRad_createAnimationForLatexReport(confidenceValue, ct, cst, slice, meanCube, mRealizations, scenProb, subIx, outpath, legendColorbar)
%%
ctDim = size(meanCube);
doseSlice = meanCube(:,:,slice);
doseSliceIx = find(doseSlice) + (slice-1)*prod(ctDim(1:2));

% check whether expDoseSliceIx is subset of subIx
assert(all(ismember(doseSliceIx, subIx)))
assert(issorted(doseSliceIx) && issorted(subIx))

[~,idxIntoSubIx] = intersect(subIx, doseSliceIx);
mRealizationsSub = mRealizations(idxIntoSubIx,:);

selectIx = doseSliceIx;
dPres = 2;
wCovMat = weightedcov(mRealizationsSub',scenProb);
mu = meanCube(selectIx);

xr = sqrt(gammaincinv(confidenceValue,numel(mu)/2)*2);


fname = 'anim';
fps = 24;
period = 5;
nFrames = period*fps;
samples = nw_getGaussianOrbitSamples(mu,wCovMat,nFrames,xr);
samplesMax = max(samples(:));
hfAnim = figure('units','normalized','outerposition',[0 0 0.5 0.5]);
%outfile = [fname '.gif'];
alpha = 0.5;
close
figure
set(gcf,'color','w');
for f=1:nFrames
    sampleCube = zeros(size(meanCube));
    sampleCube(selectIx) = samples(:,f);
    matRad_plotSliceWrapper(gca,ct,cst,1,sampleCube,3,slice,0,alpha,colorcube,jet,legendColorbar,[0.01*dPres samplesMax],[0.1 0.25 0.6 0.9 0.95 1 1.05 1.25]'*dPres);%,figXzoom,[figYzoom]);
    F(f) = getframe(gcf);
    im = frame2im(F(f));
    [imind,cm] = rgb2ind(im,256);
    
    outfile = fullfile(outpath,[fname '_' num2str(f) '.png']);
    imwrite(imind,cm,outfile,'png');
    delete(gca);
end
close
end
