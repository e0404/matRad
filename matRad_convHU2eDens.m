function eDens = matRad_convHU2eDens(HU)

% hlut table from dkfz_somatom_old

hlut = [-1024.0         0.00324;
          200.0         1.20000;
          449.0         1.20000;
         2000.0         2.49066;
         2048           2.53060;
         3071           2.53060];

eDens = interp1(hlut(:,1),hlut(:,2),HU,'linear');

if sum(isnan(eDens(:))) > 0
    error('Invalid HU range');
end