function ct = matRad_addMovement(ct, motionPeriod, numOfCtScen, amp)

if nargin < 4
    amp = [0 5 0];
end

% if isfield(ct,'dvf')
%     return
% end
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:numOfCtScen
    
    im = ct.cube{1};
    
    ct.dvf{i} = zeros([size(im), 3]);
    
    ct.dvf{i}(:,:,:,1) = amp(1) * sin(i * pi / (motionPeriod * 2));
    ct.dvf{i}(:,:,:,2) = amp(2) * sin(i * pi / (motionPeriod * 2));
    ct.dvf{i}(:,:,:,3) = amp(3) * sin(i * pi / (motionPeriod * 2));
    
    if i~=numOfCtScen 
        ct.cube{i+1} = imwarp(im, ct.dvf{i});
        
        % plotting
        subplot(121)
        imagesc(im(:,:,120))
        subplot(122)
        imagesc(ct.cube{i+1}(:,:,120))
        title(i)
        pause(.2)
        hold on
    end
    
    ct.dvf{i}(:,:,:,1) = ct.dvf{i}(:,:,:,1) * ct.resolution.x;
    ct.dvf{i}(:,:,:,2) = ct.dvf{i}(:,:,:,2) * ct.resolution.y;
    ct.dvf{i}(:,:,:,3) = ct.dvf{i}(:,:,:,3) * ct.resolution.z;
    ct.dvf{i} = permute(ct.dvf{i}, [4,1,2,3]);
    
end
close
end

