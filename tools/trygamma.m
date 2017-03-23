method='interpolation';
cube1=dose_3mm;


if exist('method','var')
    n = str2num(method(end));
    xmax=size(cube1,1)*ct.resolution.x;
    ymax=size(cube1,2)*ct.resolution.y;
    if n == 2
        Vq=zeros([size(cube1,1)*2 size(cube1,2)*2 size(cube1,3)]);
        [X,Y]=meshgrid(1 : ct.resolution.x : xmax , 1 : ct.resolution.y : ymax);
        [Xq,Yq]=meshgrid(1 : ct.resolution.x/2 : xmax , 1 : ct.resolution.y/2 : ymax);
    elseif n == 4
        Vq=zeros([size(cube1,1)*4-1 size(cube1,2)*4-1 size(cube1,3)]);
        [X,Y]=meshgrid(1 : ct.resolution.x : xmax , 1 : ct.resolution.y : ymax);
        [Xq,Yq]=meshgrid(1 : ct.resolution.x/4 : xmax , 1 : ct.resolution.y/4 : ymax);
    end
    for i=1:size(cube1,3)
        Vq(:,:,i) = interp2(X,Y,cube1(:,:,i),Xq,Yq,method(1:end-1));
    end
    
        
    if strcmp(method,'interpolation')
        Vq=zeros([size(cube1,1)*2 size(cube1,2)*2 size(cube1,3)]);
        xmax=size(cube1,1)*ct.resolution.x;
        ymax=size(cube1,2)*ct.resolution.y;
        [X,Y]=meshgrid(1 : ct.resolution.x : xmax , 1 : ct.resolution.y : ymax);
        [Xq,Yq]=meshgrid(1 : ct.resolution.x/2 : xmax , 1 : ct.resolution.y/2 : ymax);
        for i=1:size(cube1,3)
            Vq(:,:,i) = interp2(X,Y,cube1(:,:,i),Xq,Yq);
        end
    end
end
