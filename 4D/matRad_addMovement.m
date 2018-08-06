function ct = matRad_addMovement(ct, motionPeriod, numOfCtScen, amp)

if nargin < 4
    amp = [0 5 0];
end

% if there is already a dvf, no need to addMovement
if isfield(ct,'dvf') 
    return
end

figure %temp

for i = 1:numOfCtScen
    
    im = ct.cube{1};
    
    ct.dvf{i} = zeros([size(im), 3]);
    
    ct.dvf{i}(:,:,:,1) = amp(1) * sin((i - 1) * pi / (motionPeriod * 2));
    ct.dvf{i}(:,:,:,2) = amp(2) * sin((i - 1) * pi / (motionPeriod * 2));
    ct.dvf{i}(:,:,:,3) = amp(3) * sin((i - 1) * pi / (motionPeriod * 2));
    

    ct.cube{i} = imwarp(im, ct.dvf{i});
    
    % plotting %temp{
    subplot(121)
    imagesc(im(:,:,120))
    subplot(122)
    imagesc(ct.cube{i}(:,:,120))
    title(i)
    pause(.8)
    hold on  %}

    
    ct.dvf{i}(:,:,:,1) = ct.dvf{i}(:,:,:,1) * ct.resolution.x;
    ct.dvf{i}(:,:,:,2) = ct.dvf{i}(:,:,:,2) * ct.resolution.y;
    ct.dvf{i}(:,:,:,3) = ct.dvf{i}(:,:,:,3) * ct.resolution.z;
    ct.dvf{i} = permute(ct.dvf{i}, [4,1,2,3]);
    
end
close % temp
end

