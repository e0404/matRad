function [finalWeight, sigmaBeamlet, posX, posY, numOfSub] = ...
                        matRad_calcWeights(sigmaTot, method, N, sigmaSub)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates weights for a fine sampling pencil beam
% algorithm
%
% call
%   [finalWeight, sigmaBeamlet, posX, posY, numOfSub] = ...
%                       matRad_calcWeights(sigmaTot, method, N, sigmaSub)
%
% input
%   sigmaTot:       the standard deviation of the lateral spread of the pencil
%                   beam
%   method:         method of weight calculation
%                   'russo' for integral calculation according to [1]
%                   'fitCircle'   for using square grid weights derived from a fit
%                   'fitSquare'   for using circular grid weights derived from a fit
%   N:              if method == russo:
%                   number of subsample beams shells. Means we have a
%                   grid of NxN sub beams representing the total beam
%                   if method == fitCircle or fitSquare:
%                   number of subsample beams shells. n = 2 means we have two 
%                   lines of points around the central ray. The number of 
%                   sub-beams will be:
%                   #sb = (2*n +1)^2 for the square;
%                   #sb = (2^n -1)*6 +1 for the circle
%   sigmaSub:       is the sigma of the gaussian of the sub-beams (only
%                   needed for russo method
%
% output
%   finalWeight:    is the array of the weights of the sub-pencil beams. It
%                   runs over the same index as posx and posy
%   
%   posx & posy:    are the positions of the sub-pencil beams, returned as
%                   meshgrid-matrices if method is 'square' and as vectors 
%                   if method is 'circle'
%   numOfSub:       number of sub-pencil beams
%
% References
%   [1] https://iopscience.iop.org/article/10.1088/0031-9155/61/1/183
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2020 the matRad development team.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~strcmp(method, 'russo') && ~strcmp(method, 'fitCircle') && ~strcmp(method, 'fitSquare')
    
    error('method not supported');

elseif strcmp(method, 'russo')
    % splitting into N^2 sub beams accoring to Russo et al (2016)
    
    sigmaHead = sqrt(sigmaTot^2 - sigmaSub^2);

    gauss = @(sigma, x, y, muX, muY) 1 / (2 * pi * sigma^2) .* exp(-((x + muX).^2 + (y + muY).^2) / (2 * sigma^2)); 

    R = 3.5 * sigmaHead;
    dR = 2 * R / N;

    counter = 1;
    for iX = -(N - 1) / 2 : (N - 1) / 2
        for iY = -(N - 1) / 2 : (N - 1) / 2
            finalWeight(counter) = integral2(@(x,y) gauss(sigmaHead, x, y, 0, 0), ...
                            (iX - 0.5) * dR, (iX + 0.5) * dR, ...
                            (iY - 0.5) * dR, (iY + 0.5) * dR);
            posX(counter) = iX * dR;
            posY(counter) = iY * dR;
            sigmaBeamlet(counter) = sigmaSub;

            counter = counter + 1;
        end
    end

    finalWeight = finalWeight';
    posX = posX';
    posY = posY';
    sigmaBeamlet = sigmaBeamlet';


    numOfSub = N * N;
    
elseif strcmp(method, 'fitCircle') || strcmp(method, 'fitSquare')
    % number of sub beams will be (2*n +1)^2 for the square;
    %                             (2^n -1)*6 +1 for the circle
    if N~=2 && N~=3 && N~=8
        error('number of shells N not supported');
    end

    % This parameters come from simulations done previously
    % see "Research on the dosimetric accuracy of pencil beam fine sampling
    % for radiation proton dose calculation" by Giuseppe Pezzano (2018)
    if N == 2
        if strcmp(method,'fitCircle')
            sigmaBeamlet = 0.8237 .* sigmaTot;
            radius    = 0.6212 .* sigmaTot;
            X1(1,:)     = 0.3866 .* sigmaTot.^2;
            X1(2,:)     = 0.6225 .* sigmaTot;
        elseif strcmp(method,'fitSquare')
            sigmaBeamlet = 0.8409 .* sigmaTot;
            radius    = 0.5519 .* sigmaTot;
            X1(1,:)     = 0.3099 .* sigmaTot.^2;
            X1(2,:)     = 0.5556 .* sigmaTot;
        end
    elseif N == 3
        if strcmp(method,'fitCircle')
            sigmaBeamlet = 0.7605 .* sigmaTot;
            radius    = 0.5000 .* sigmaTot;
            X1(1,:)     = 0.3006 .* sigmaTot.^2 - 1.3005 .* sigmaTot + 7.3097;
            X1(2,:)     = 0.6646 .* sigmaTot - 0.0044;
        elseif strcmp(method,'fitSquare')
            sigmaBeamlet = 0.8409 .* sigmaTot;
            radius    = 0.5391 .* sigmaTot + 0.0856;
            X1(1,:)     = 0.3245 .* sigmaTot.^2 + 0.0001 .* sigmaTot - 0.0004;
            X1(2,:)     = 0.6290 .* sigmaTot - 0.0403;
        end
    elseif N == 8 && strcmp(method,'fitCircle')
        sigmaBeamlet = 0.5 .* sigmaTot;
        radius    = 0.25 .* sigmaTot;
        X1(1,:)     = 0.0334 .* sigmaTot.^2 - 4.1061e-06 .* sigmaTot + 1.5047e-06;
        X1(2,:)     = 0.6 .* sigmaTot + 3.3151e-06;
    else
        error('number of shells N not supported');
    end

    % setting positions of sub-beams
    if strcmp(method,'fitSquare')
        numOfSub = (2*N +1)^2;
        points   = linspace(-radius*N,radius*N,2*N +1);
        posX     = points'*ones(1,2*N +1);
        posY     = posX';
    else
        dim = size(radius,2);
        numOfSub = (2^N -1)*6 +1;
        ang  = zeros(1,1);
        posX = zeros(1,dim);
        posY = zeros(1,dim);
        radiusShell = zeros(1,dim);
        for i = 1:N
            subsInShell = 6 * 2^(i-1);
            % this takes the sub-beams index in one shell
            ang         = cat(2, ang, pi .* linspace(0,2-2/subsInShell, subsInShell));
            radiusShell = cat(1, radiusShell, ones(subsInShell,1)*(i.*radius));
        end
        posX = cat(1, posX, bsxfun(@times,cos(ang(2:end))',radiusShell(2:end,:)));
        posY = cat(1, posY, bsxfun(@times,sin(ang(2:end))',radiusShell(2:end,:)));
    end

    % compute weights at positions
    sig  = ones(size(posX,1),1)*X1(2,:);
    normSig = ones(size(posX,1),1)*X1(1,:);

    finalWeight = normSig .* (2.*pi.*sig.^2).^(-1) .* exp(-(posX.^2+posY.^2)./(2.*(sig.^2)));
    finalWeight = reshape(finalWeight, numel(finalWeight), 1);
    sigmaBeamlet = repmat(reshape(sigmaBeamlet, numel(sigmaBeamlet), 1), numel(finalWeight),1);
    posX =reshape(posX, numel(posX), 1);
    posY =reshape(posY, numel(posY), 1);
    
end
