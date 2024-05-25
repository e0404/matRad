function [gammaPassRate,gammaValues] =...
    matRad_gammaPassRate2D(doseMes,doseCalc, grids, DTA , doseDiff, frac )

% matRad_gammaPassRate2D compares two dose distributions in terms of gamma
% analysis.
%
% call
%   [gammaPassRate,gammaValues] =...
%    matRad_gammaPassRate2D(doseMes,doseCalc, grids, DTA , doseDiff, frac )
%
% input
%   doseMes: Measured dose (reference), 3D dose distribution grid
%   doseCalc: Calculated dose (to evaluate, -"-)
%   grids: struct
%       grids.x: 3D meshgrid with x coords
%       grids.y: 3D meshgrid with y coords 
%       grids.z: 3D meshgrid with y coords 
%   DTA: Distance To Agreement in mm (criterion parameter)
%   doseDiff: doseDiff in % (criterion parameter)
%   frac: fraction of random dose points evaluated
%         (this is used to safe time while get a statistical result)
%
% output
%   3D grid with gamma factors
%
% reference
%   see D. A. Low, W. B. Harms, S. Mutic, and J. A. Purdy:
%       „A technique for the quantitativeevaluation of dose distributions“,
%       Medical Physics, vol. 25, no. 5, pp. 656–661, 1998.
%       DOI:https://doi.org/10.1118/1.598


% define and explain quantities...
% d: distance to agreement [mm]
d = DTA;

% D: dose distance as fraction from measured dose
D = doseDiff/100;

% preallocation for grid of all gamma values
gammaValues = zeros(size(doseCalc));

% calculate gamma accordance for all measurement points at positions rm
% with array indey rmIndex:

% initialize waitbar
w = waitbar(0,'validating individual points');

% create random fraction
indices = 1:numel(doseMes);
randFrac = rand(size(indices)) < frac;
usedIndices = indices(randFrac);
run = 1;

%run for loop
for rmIndex = usedIndices
    % rm: r measurement vector
    rm = [grids.x(rmIndex), grids.y(rmIndex), grids.z(rmIndex)]; 
    
    Gamma = sqrt(((grids.x-rm(1)).^2 + (grids.y-rm(2)).^2 +...
        (grids.z-rm(3)).^2)/d^2 + (doseMes(rmIndex) - doseCalc).^2/D^2); 
              
    gamma = min(Gamma,[],'all');
    
    gammaValues(rmIndex) = gamma;
    
    waitbar(run/numel(usedIndices),w);
    run = run+1;
end

gammaPassRate = nnz(gammaValues(usedIndices)<1)/numel(usedIndices);

close(w)

end