function [sigmaNoRashi, sigmaWithRashi] = helper_rashiLateralBroadening(engineName, propDoseCalc, testDataFile)
% helper_rashiLateralBroadening  Measures lateral beam broadening due to a
% range shifter for a given particle pencil-beam dose engine.
%
% Runs a forward dose calculation on the protons test dataset (first beam
% only) once without and once with a thick range shifter on all rays, and
% returns the lateral beam width (second moment of the profile) at a
% shallow depth slice for both cases. The range shifter is placed with a
% large air gap so that its scattering contribution clearly dominates the
% lateral width (see issue #923).
%
% call
%   [sigmaNoRashi, sigmaWithRashi] = helper_rashiLateralBroadening(engineName)
%   [sigmaNoRashi, sigmaWithRashi] = helper_rashiLateralBroadening(engineName, propDoseCalc)
%   [sigmaNoRashi, sigmaWithRashi] = helper_rashiLateralBroadening(engineName, propDoseCalc, testDataFile)
%
% inputs
%   engineName   shortName of the dose engine to use (e.g. 'HongPB')
%   propDoseCalc (optional) struct with additional pln.propDoseCalc fields
%   testDataFile (optional) name of the test dataset to use
%                (default 'protons_testData.mat')
%
% outputs
%   sigmaNoRashi   lateral profile sigma [mm] without range shifter
%   sigmaWithRashi lateral profile sigma [mm] with range shifter

if nargin < 3
    testDataFile = 'protons_testData.mat';
end

testData = load(testDataFile);

pln = testData.pln;
pln.propDoseCalc.engine = engineName;
if nargin > 1 && ~isempty(propDoseCalc)
    fn = fieldnames(propDoseCalc);
    for i = 1:numel(fn)
        pln.propDoseCalc.(fn{i}) = propDoseCalc.(fn{i});
    end
end

% use only the first beam (gantry & couch angle 0, i.e. beam travels along +y)
stf = testData.stf(1);

noRashi = struct('ID', 0, 'eqThickness', 0, 'sourceRashiDistance', 0);
% thick range shifter close to the source (large air gap) for a strong
% scattering contribution
rashi   = struct('ID', 1, 'eqThickness', 40, 'sourceRashiDistance', 400);

stfNoRashi = stf;
stfWithRashi = stf;
for r = 1:numel(stf.ray)
    nEnergies = numel(stf.ray(r).energy);
    stfNoRashi.ray(r).rangeShifter   = repmat(noRashi, 1, nEnergies);
    stfWithRashi.ray(r).rangeShifter = repmat(rashi, 1, nEnergies);
end

w = ones(stf.totalNumOfBixels, 1);
resultNoRashi   = matRad_calcDoseForward(testData.ct, testData.cst, stfNoRashi, pln, w);
resultWithRashi = matRad_calcDoseForward(testData.ct, testData.cst, stfWithRashi, pln, w);

% evaluate the lateral (x) second moment at a shallow slice along the beam
% (+y) direction, where the depth-dependent MCS contribution is small
sliceIx = 2;
sigmaNoRashi   = lateralSigma(squeeze(resultNoRashi.physicalDose(sliceIx, :, :)), testData.ct.resolution.x);
sigmaWithRashi = lateralSigma(squeeze(resultWithRashi.physicalDose(sliceIx, :, :)), testData.ct.resolution.x);
end

function s = lateralSigma(slice2D, res)
% profile-weighted standard deviation along x of a 2D (x,z) dose slice
profile = sum(slice2D, 2);
x = (1:numel(profile))' * res;
mu = sum(profile .* x) / sum(profile);
s = sqrt(sum(profile .* (x - mu).^2) / sum(profile));
end
