function Untitled2(sigma_ray,sigma_sub,radius,n,X1,method)

clear subGauss

if ~exist('method','var')
    method = 'square';
end

if ~strcmp(method,'square') && ~strcmp(method,'circle')
    error('method not supported');
end

% setting positions of sub-beams
if strcmp(method,'square')
    numOfSub = (2*n +1)^2;
    points = linspace(-radius*(sqrt(numOfSub)-1)/2,radius*(sqrt(numOfSub)-1)/2,sqrt(numOfSub));
    posx = points'*ones(1,sqrt(numOfSub));
    posy = posx';
else
    numOfSub = (2^n -1)*6 +1;
    ang = zeros(1,1);
    posx = zeros(1,1);
    posy = zeros(1,1);
    radiusShell = zeros(1,1);
    for i=1:n
        SubsInShell = (2^i -1)*6 +1 - ((2^(i-1) -1)*6 +1 );
        % this takes the sub-beams index in one shell
        ang = cat(2, ang, pi .* linspace(0,2-2/SubsInShell, SubsInShell));
        radiusShell = cat(2, radiusShell, i.*radius.*ones(1, SubsInShell));
    end
    posx = cat(2, posx, posx(1) + radiusShell(2:end).*cos(ang(2:end)));
    posy = cat(2, posy, posy(1) + radiusShell(2:end).*sin(ang(2:end)));
end

x0 = -3*sigma_ray:sigma_ray/70:3*sigma_ray;
y0 = x0;


gaussian2 = @(x, y, mux, muy ,sig) (2*pi*sig^2)^(-1) .* exp(-(x-mux).^2/(2*(sig^2)))' * exp(-(y-muy).^2/(2*(sig^2)));
f1test = gaussian2(x0,y0,0,0,sigma_ray);
f2test = @(X) 0;
for i=1:numOfSub
    f2test = @(X) f2test(X) + X(1) .* ...
        gaussian2(posx(i),posy(i),0,0,X(2)).*gaussian2(x0,y0,posx(i),posy(i),sigma_sub);
    subGauss(:,:,i) = X1(1) .* ...
        gaussian2(posx(i),posy(i),0,0,X1(2)).*gaussian2(x0,x0,posx(i),posy(i),sigma_sub);
end

maxibon = max(max(abs((f2test(X1)-f1test)./max(max(f1test)).*100)));

figure
% subplot(2,1,1)
%hold off
%scatter3(0,0,0)
hold
if n==2
    if method=='square'
        colorvec = ['y','y','y','y','y','y','g','c','g','y','y','c','b','c','y','y','g','c','g','y','y','y','y','y','y'];
        falpha = [1, 1, 1, 1, 1, 1, .8, .4, .8, 1, 1, .4, .2, .4, 1, 1, .8, .4, .8, 1, 1, 1, 1, 1, 1];
    else
        colorvec = ['b','g','c','g','c','g','c','y','y','y','y','y','y','y','y','y','y','y','y'];
        falpha = [.2, .4, .4, .4, .4, .4, .4, .8, .8, .8, .8, .8, .8, .8, .8, .8, .8, .8, .8];
    end
else if n==3
        if method=='square'
            colorvec = ['m','m','m','m','m','m','m','m','y','y','y','y','y','m','m','y','g','c','g','y','m','m','y','c','b','c','y','m','m','y','g','c','g','y','m','m','y','y','y','y','y','m','m','m','m','m','m','m','m'];
            falpha = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, .8, .4, .8, 1, 1, 1, 1, .4, .2, .4, 1, 1, 1, 1, .8, .4, .8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
        else
            colorvec = ['b','g','c','g','c','g','c','y','y','y','y','y','y','y','y','y','y','y','y','m','m','m','m','m','m','m','m','m','m','m','m','m','m','m','m','m','m','m','m','m','m','m','m'];
            falpha = [.2, .4, .4, .4, .4, .4, .4, .8, .8, .8, .8, .8, .8, .8, .8, .8, .8, .8, .8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
        end
    end
end

for i=1:numOfSub
    surf(subGauss(:,:,i),'LineWidth',1,'EdgeAlpha',.00001,'FaceAlpha',falpha(i),'FaceColor',colorvec(i))
end
title(strcat('r = ', num2str(radius), '       \sigma_s = ', num2str(sigma_sub), ...
    '         \sigma_t = ', num2str(sigma_ray)))
view([66 15])
%axis([100 350 100 350])

figure
% subplot(2,1,2)
surf(((f2test(X1)-f1test)./max(max(f1test)).*100),'LineWidth',.1,'EdgeAlpha',.3)
title(strcat('max error = ', num2str(maxibon),' %','        (max-min) diff = ',...
    num2str((max(max((f2test(X1)-f1test)./max(max(f1test)).*100))- min(min((f2test(X1)-f1test)./max(max(f1test)).*100)))/max(max(abs((f2test(X1)-f1test)./max(max(f1test)).*100)) )*100 ),...
    ' %       aver = ', num2str(mean(mean( (f2test(X1)-f1test)./max(max(f1test)).*100 /max(max(abs((f2test(X1)-f1test)./max(max(f1test)).*100)) )))) ))
view(3)