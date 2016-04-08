function cst = matRad_addTargetRing(cst,ct,ringsize)

% find all target voxels from cst cell array
V = [];
for i=1:size(cst,1)
    if isequal(cst{i,3},'TARGET') && ~isempty(cst{i,6})
        V = [V;vertcat(cst{i,4}{:})];
    end
end

% Remove double voxels
V = unique(V);

% generate voi cube for targets
voiTarget    = zeros(ct.cubeDim);
voiTarget(V) = 1;

myMargin.x = ringsize;
myMargin.y = ringsize;
myMargin.z = ringsize;
voiTarget  = matRad_addMargin(voiTarget,cst,ct.resolution,myMargin,true);
VwithMargin = find(voiTarget>0);

% add TARGET RING structure to cst
cstNewLine = size(cst,1) + 1 ;
targetLine = 2;

cst{cstNewLine,1} = cst{cstNewLine - 1,1} + 1;
cst{cstNewLine,2} = 'TargetRing';
cst{cstNewLine,3} = 'TARGET';

cst{cstNewLine,4}{1} = VwithMargin(~ismember(VwithMargin,V));

cst{cstNewLine,5} = cst{targetLine,5};    
cst{cstNewLine,6}.type = 'min DCH objective';
cst{cstNewLine,6}.penalty = 100;
cst{cstNewLine,6}.dose = 57.5;
cst{cstNewLine,6}.EUD = NaN;
cst{cstNewLine,6}.volume = 0.9;
cst{cstNewLine,6}.coverage = 0.9;
cst{cstNewLine,6}.robustness = 'coverage';

end