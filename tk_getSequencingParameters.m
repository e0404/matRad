function [x_l, x_r, shapeIx] = tk_getSequencingParameters(seqResult,pln)

% initializing variables
x_l = [];
x_r = [];
shapeIx = [];
shapeCounter = 0;


% loop over all beams
for i=1:pln.numOfBeams
    
%   loop over all shapes
    for j=1:seqResult.beam(i).numOfShapes
        shapeCounter = shapeCounter +1; % get index of current shape
        tempShape = seqResult.beam(i).shapes(:,:,j);
        [Ix_l, Ix_r] = tk_getLeafPos(tempShape);
        
        x_l = [x_l; Ix_l];
        x_r = [x_r; Ix_r];
        shapeIx = [shapeIx; shapeCounter*ones(size(Ix_l,1),1)];
        
        
        
    end
end
    
end