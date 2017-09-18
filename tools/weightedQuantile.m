function [wQ] = weightedQuantile(values, percentiles, weight, sortedB, extraPolB)

  % values vector
  % quantiles what you want to know
  % weight vector of same length as values
  % sortedB optional
  % extrapolation alloud
  if ~exist('sortedB', 'var') || isempty(sortedB)
    sortedB = false;
  end
  if ~exist('weight', 'var') || isempty(weight)
    weight = ones(size(values));
  end
  if ~exist('extraPolB', 'var') || isempty(extraPolB)
    extraPolB = false;
  end
  if ~(percentiles(:) >= 0 & percentiles(:) <= 1)
    error('Quantiles shold be [0, 1]');
  end

  if sortedB == false
    [values, sortIx] = sort(values);
    weight = weight(sortIx);
  end
  wQtemp = cumsum(weight) - 0.5 * weight;
  wQtemp = wQtemp ./ sum(weight);

  wQ = NaN * ones(size(values,1), 2);
  
  if extraPolB
      [x,ia] = unique(wQtemp);
     
      V = values(ia);
      F = griddedInterpolant(x, V);
      wQ = F(percentiles);
  else

  end
  
  wQ = wQ';

end % eof
