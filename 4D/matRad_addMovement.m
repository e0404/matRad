function ct = matRad_addMovement(ct, motionPeriod, amp)

if nargin < 3
    amp = [5 0 0];
end

for i = 1:ct.numOfCtScen
    
    im = ct.cube{1};
    
    ct.dvf{i} = zeros([size(im), 3]);
    
    ct.dvf{i}(:,:,:,1) = amp(1) * sin(i * pi / (motionPeriod * 2));
    ct.dvf{i}(:,:,:,2) = amp(2) * sin(i * pi / (motionPeriod * 2));
    ct.dvf{i}(:,:,:,3) = amp(3) * sin(i * pi / (motionPeriod * 2));
    
    if i~=ct.numOfCtScen
        ct.cube{i+1} = imwarp(im, ct.dvf{i});
    end
    
end
