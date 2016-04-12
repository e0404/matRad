function c = matRad_constFunc(d_i,constraint,d_ref)


numOfVoxels = numel(d_i);

if isequal(constraint.type, 'max dose constraint')

    epsilon = 1e-3;
    d_i_max = max(d_i);

    c = d_i_max + epsilon * log( sum(exp((d_i - d_i_max)/epsilon)) );

elseif isequal(constraint.type, 'min dose constraint')

    epsilon = 1e-3;
    d_i_min = min(d_i);

    c = d_i_min - epsilon * log( sum(exp((d_i_min - d_i)/epsilon)) );

elseif isequal(constraint.type, 'min mean dose constraint') || ...
       isequal(constraint.type, 'max mean dose constraint') || ...
       isequal(constraint.type, 'min max mean dose constraint')

    c = mean(d_i);

elseif isequal(constraint.type, 'min EUD constraint') || ...
       isequal(constraint.type, 'max EUD constraint') || ...
       isequal(constraint.type, 'min max EUD constraint')

    exponent = constraint.EUD;

    c = mean(d_i.^exponent)^(1/exponent);

elseif isequal(constraint.type, 'exact DVH constraint') || ...
       isequal(constraint.type, 'max DVH constraint') || ... 
       isequal(constraint.type, 'min DVH constraint')

    c = sum(d_i >= d_ref)/numOfVoxels ;

    % alternative constraint calculation 3/4 %
    % % get reference Volume
    % refVol = cst{j,6}(k).volume/100;
    % 
    % % calc deviation
    % deviation = d_i - d_ref;
    % 
    % % calc d_ref2: V(d_ref2) = refVol
    % d_ref2 = matRad_calcInversDVH(refVol,d_i);
    % 
    % % apply lower and upper dose limits
    % if isequal(cst{j,6}(k).type, 'max DVH constraint')
    %    deviation(d_i < d_ref | d_i > d_ref2) = 0;
    % elseif isequal(cst{j,6}(k).type, 'min DVH constraint')
    %    deviation(d_i > d_ref | d_i < d_ref2) = 0;
    % end
    % 
    % %c = sum(deviation);                              % linear deviation
    % %c = deviation'*deviation;                        % square devioation
    % c = (1/size(cst{j,4},1))*(deviation'*deviation); % square deviation with normalization
    % %c = (deviation).^2'*(deviation).^2;               % squared square devioation
    % alternative constraint calculation 3/4 %

end 


end