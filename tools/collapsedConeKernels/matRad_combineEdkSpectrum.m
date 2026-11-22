function [poly, spectrum] = matRad_combineEdkSpectrum(mono, spectrumFile)
% matRad_combineEdkSpectrum combines monoenergetic energy deposition kernels
% into a polyenergetic kernel weighted by the TERMA fraction of each energy
%
% call:
%   [poly, spectrum] = matRad_combineEdkSpectrum(mono, spectrumFile)
%
% input:
%   mono:           struct array of monoenergetic kernels as returned by
%                   matRad_readEdkListing, sharing the same geometry
%   spectrumFile:   EGSnrc spectrum file (e.g. mohan6.spectrum)
%
% output:
%   poly:           polyenergetic kernel with the same fields as the
%                   monoenergetic kernels, normalized per unit TERMA
%   spectrum:       struct with the spectral bin centers [MeV], fluence,
%                   water attenuation coefficients [cm^2/g] and the
%                   normalized TERMA weights of each bin
%
% References
%   [1] Mohan, Chui, Lidofsky, Med. Phys. 12, 592 (1985)
%   [2] Hubbell, Seltzer, NISTIR 5632 (1995)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matRad_cfg = MatRad_Config.instance();

[~, order] = sort([mono.energy]);
mono = mono(order);
monoEnergies = [mono.energy];

for i = 2:numel(mono)
    if ~isequal(mono(i).angleEdges, mono(1).angleEdges) || ~isequal(mono(i).radiusEdges, mono(1).radiusEdges)
        matRad_cfg.dispError('All monoenergetic kernels must share the same geometry.');
    end
end

spectrum = matRad_readEgsSpectrum(spectrumFile);
spectrum.muRho = matRad_waterAttenuation(spectrum.energy);
spectrum.termaWeight = spectrum.fluence .* spectrum.energy .* spectrum.muRho;
spectrum.termaWeight = spectrum.termaWeight / sum(spectrum.termaWeight);

monoWeight = zeros(1, numel(mono));
eClamped = min(max(spectrum.energy, monoEnergies(1)), monoEnergies(end));
for j = 1:numel(eClamped)
    iHigh = find(monoEnergies >= eClamped(j), 1, 'first');
    if monoEnergies(iHigh) == eClamped(j)
        monoWeight(iHigh) = monoWeight(iHigh) + spectrum.termaWeight(j);
    else
        iLow = iHigh - 1;
        f = (eClamped(j) - monoEnergies(iLow)) / (monoEnergies(iHigh) - monoEnergies(iLow));
        monoWeight(iLow) = monoWeight(iLow) + (1 - f) * spectrum.termaWeight(j);
        monoWeight(iHigh) = monoWeight(iHigh) + f * spectrum.termaWeight(j);
    end
end
spectrum.monoEnergies = monoEnergies;
spectrum.monoWeight = monoWeight;

poly = struct();
poly.energy = sum(monoWeight .* monoEnergies);
poly.numHistories = sum([mono.numHistories]);
poly.angleEdges = mono(1).angleEdges;
poly.radiusEdges = mono(1).radiusEdges;

comps = {'total', 'primary', 'scatter'};
for c = 1:numel(comps)
    val = zeros(size(mono(1).total));
    varAbs = zeros(size(mono(1).total));
    for i = 1:numel(mono)
        v = mono(i).(comps{c});
        u = mono(i).([comps{c} 'RelUnc']) .* v;
        u(~isfinite(u)) = 0;
        val = val + monoWeight(i) * v;
        varAbs = varAbs + (monoWeight(i) * u).^2;
    end
    relUnc = inf(size(val));
    relUnc(val > 0) = sqrt(varAbs(val > 0)) ./ val(val > 0);
    poly.(comps{c}) = val;
    poly.([comps{c} 'RelUnc']) = relUnc;
end

end

function spectrum = matRad_readEgsSpectrum(fileName)
fid = fopen(fileName, 'r');
if fid < 0
    matRad_cfg = MatRad_Config.instance();
    matRad_cfg.dispError('Could not open spectrum file %s.', fileName);
end
cleaner = onCleanup(@() fclose(fid));
spectrum.title = strtrim(fgetl(fid));
header = sscanf(strrep(fgetl(fid), ',', ' '), '%f');
data = fscanf(fid, '%f , %f', [2 Inf]);

numBins = header(1);
eMin = header(2);
specMode = header(3);
upperEdge = data(1, 1:numBins);
values = data(2, 1:numBins);
lowerEdge = [eMin upperEdge(1:end - 1)];
width = upperEdge - lowerEdge;

spectrum.energy = 0.5 * (lowerEdge + upperEdge);
if specMode == 1
    spectrum.fluence = values .* width;
else
    spectrum.fluence = values;
end
spectrum.fluence = spectrum.fluence / sum(spectrum.fluence);
end

function muRho = matRad_waterAttenuation(energy)
tab = [0.10 0.1707
       0.15 0.1505
       0.20 0.1370
       0.30 0.1186
       0.40 0.1061
       0.50 0.09687
       0.60 0.08956
       0.80 0.07865
       1.00 0.07072
       1.25 0.06323
       1.50 0.05754
       2.00 0.04942
       3.00 0.03969
       4.00 0.03403
       5.00 0.03031
       6.00 0.02770
       8.00 0.02429
       10.0 0.02219];
muRho = exp(interp1(log(tab(:, 1)), log(tab(:, 2)), log(energy), 'linear', 'extrap'));
end
