function S = wMean(X,w)
if exist('w','var') || ~isempty(w)
    if isvector(X) && isvector(w)
        S = reshape(w,1,[]) * reshape(X,[],1) / sum(w);
    else
        % row-wise
        S = reshape(w,1,[]) * X ./ sum(w);        
    end
    
else
    S = mean(X);
end
end

