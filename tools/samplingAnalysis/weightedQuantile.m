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
      x = single(wQtemp);
      V = values;
      %samplePoints = {x, single(1:size(V,1))};
      %queryPoints = {single(percentiles), single(1:size(V,1))};
      F = griddedInterpolant(x, V);
      wQ = F(percentiles);
  else

  end
  
  wQ = wQ';
  
  % give the ranger of the most outer quantiles
  %ix = NaN;
  %ix = getRange(values, [wQ(wQ == min(wQ(:))), wQ(wQ == max(wQ(:)))]);

end % eof
