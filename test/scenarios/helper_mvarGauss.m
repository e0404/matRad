function p = helper_mvarGauss(model)
%helper_mvarGauss Computes multivariate Gaussian probability for scenario
%model
%   Can be used to test correct scenario probabilities. Also considers used
%   ct phases for probability computation (uncorrelated)
    Sigma = diag([model.shiftSD,model.rangeAbsSD,model.rangeRelSD./100].^2);
    d = size(Sigma,1);
    [cs,~] = chol(Sigma);
    
    % Compute from Gaussian errors
    p = (2*pi)^(-d/2) * exp(-0.5*sum((model.scenForProb(:,2:end)/cs).^2, 2)) / prod(diag(cs));

    % Now multiplay with the phase probability
    tmpPhaseProb = arrayfun(@(phase) model.ctScenProb(find(model.ctScenProb(:,1) == phase),2),model.scenForProb(:,1));
    p = p .* tmpPhaseProb;
end