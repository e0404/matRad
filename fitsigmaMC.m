load protons_generic_MCsquare.mat
sigma = [];
for i = 1:numel(machine.data)
    newDepths = linspace(0,machine.data(i).depths(end),numel(machine.data(i).depths) * 100);
    newDose   = interp1(machine.data(i).depths, machine.data(i).Z, newDepths, 'spline');
     
    machine.data(i).depths = newDepths;
    machine.data(i).Z      = newDose;
    
    %interpolate range at 80% dose after peak.
    [maxV, maxI] = max(machine.data(i).Z);
    [~, r80ind] = min(abs(machine.data(i).Z(maxI:end) - 0.8 * maxV));
 
    %find FWHM w50 of bragg peak
    [~, d50rInd] = min(abs(machine.data(i).Z(maxI:end) - 0.5 * maxV));
    d50rInd = d50rInd - 1;
    d50_r = interp1(machine.data(i).Z(maxI + d50rInd - 1:maxI + d50rInd + 1), ...
                            machine.data(i).depths(maxI + d50rInd - 1:maxI + d50rInd + 1), 0.5 * maxV);
                             
    array1 = machine.data(i).depths(maxI-round(d50rInd):maxI+d50rInd);
    array2 = machine.data(i).Z(maxI-round(d50rInd):maxI+d50rInd);  

    gauss = @(x, a, sigma, b) a * exp(- (x - b).^2 / (2*sigma^2));

    funcs.objective = @(p) sum((gauss(array1, p(1), p(2),p(3)) - array2).^2);

    funcs.gradient = @(p) [ sum(2 * (p(1) * exp(- (array1 - p(3)).^2 / (2 * p(2)^2)) - array2) ...
                                            .* exp(- (array1 - p(3)).^2 / (2 * p(2)^2)));
                            sum(2 * (p(1) * exp(- (array1 - p(3)).^2 / (2 * p(2)^2)) - array2) ...
                                            .* p(1) .* exp(- (array1 - p(3)).^2 / (2 * p(2)^2)) .* (array1 - p(3)).^2 / p(2)^3)
                            sum(2 * (p(1) * exp(- (array1 - p(3)).^2 / (2 * p(2)^2)) - array2) ...
                                            .* p(1) .* exp(- (array1 - p(3)).^2 / (2 * p(2)^2)) .* 2 .* (array1 - p(3)) / (2 * p(2)^2))];


    options.lb = [0,  0, 0];
    options.ub = [ Inf, Inf, Inf];
    options.ipopt.limited_memory_update_type = 'bfgs';
            options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.print_level = 1;

    start = [maxV, 2 * d50_r, machine.data(i).depths(maxI)];
   
    [fitResult, ~] = ipopt (start, funcs, options);

    sigma = [sigma, fitResult(2)];
end