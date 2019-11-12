clear
load fitPointsGeneric.mat
sigmasG = sigEneSpr(1,:);
eneG    = sigEneSpr(2,:);
bestG   = sigEneSpr(3,:);
load fitPointsTopas.mat
sigmasT = sigEneSpr(1,:);
bestT   = sigEneSpr(3,:);
eneT    = sigEneSpr(2,:);
load fitPointsMCs.mat
sigmasM = rere(1,:);
bestM   = rere(3,:);
eneM    = rere(2,:);
load fitPointsBDLmat.mat
sigmasB = spread(1,:);
bestB   = spread(3,:);
eneB    = spread(2,:);

comSig  = [sigmasG, sigmasT, sigmasM, sigmasB]';
comEne  = [eneG, eneT, eneM, eneB]';
comBest = [bestG, bestT, bestM, bestB]';

% fitFunction = @(E, sigma, A1, A2, k1, k2, k3, b1, b2) ...
%                     A1 .* (1 - exp(-k1 .* (E - b1))) .* (1 - 1 ./ ( 1 + exp(-k3 .* (sigma-b2)))) + ...
%                     A2 .*      exp(-k2 .* (E - b1))  .*      1 ./ ( 1 + exp(-k3 .* (sigma-b2)));
%                 
% funcs.objective = @(p) sum((fitFunction(comEne, comSig, p(1), p(2), p(3), p(4), p(5), p(6), p(7)) - comBest).^2);
% 
% funcs.gradient = @(p) [ 2 * sum((fitFunction(comEne, comSig, p(1), p(2), p(3), p(4), p(5), p(6), p(7)) - comBest) .* ...
%                             (1 - exp(-p(3) .* (comEne - p(6)))) .* (1 - 1 ./(1 + exp(-p(5) .* (comSig - p(7))))));
%                         2 * sum((fitFunction(comEne, comSig, p(1), p(2), p(3), p(4), p(5), p(6), p(7)) - comBest) .* ...
%                             exp(-p(4) .* (comEne - p(6))) ./ (1 + exp(-p(5) .* (comSig - p(7)))));
%                         2 * sum((fitFunction(comEne, comSig, p(1), p(2), p(3), p(4), p(5), p(6), p(7)) - comBest) .* ...
%                             p(1) .* (comEne - p(6)) .* exp(-p(3) .* (comEne - p(6))) .* (1 - (1 ./ (1 + exp(-p(5) .* (comSig - p(7)))))));
%                         2 * sum((fitFunction(comEne, comSig, p(1), p(2), p(3), p(4), p(5), p(6), p(7)) - comBest) .* ...
%                             p(2) .* -(comEne - p(6)) .* exp(-p(4) .* (comEne - p(6))) ./ (1 + exp(-p(5) .* (comSig - p(7)))));
%                         2 * sum((fitFunction(comEne, comSig, p(1), p(2), p(3), p(4), p(5), p(6), p(7)) - comBest) .* ...
%                             p(1) .* (1 - exp(-p(3) .* (comEne - p(6)))) .* -(comSig - p(7)) .* exp(-p(5) .* (comSig - p(7))) ./ ((1 + exp(-p(5) .* (comSig - p(7)))).^2) + ...
%                             p(2) .* exp(-p(4) .* (comEne - p(6))) .* (comSig - p(7)) .* exp(-p(5) .* (comSig - p(7))) ./ ((1 + exp(-p(5) .* (comSig - p(7)))).^2));
%                         2 * sum((fitFunction(comEne, comSig, p(1), p(2), p(3), p(4), p(5), p(6), p(7)) - comBest) .* ...
%                             p(1) .* -p(3) .* exp(-p(3) .* (comEne - p(6))) .* (1 - 1 ./(1 + exp(-p(5) .* (comSig - p(7))))) + ...
%                             p(2) .*  p(4) .* exp(-p(4) .* (comEne - p(6))) ./ (1 + exp(-p(5) .* (comSig - p(7)))));
%                         2 * sum((fitFunction(comEne, comSig, p(1), p(2), p(3), p(4), p(5), p(6), p(7)) - comBest) .* ...
%                             p(1) .* (1 - exp(-p(3) .* (comEne - p(6)))) .* p(5) .* exp(-p(5) .* (comSig - p(7))) ./ ((1 + exp(-p(5) .* (comSig - p(7)))).^2) + ...
%                             p(2) .* exp(-p(4) .* (comEne - p(6))) .* -p(5) .* exp(-p(5) .* (comSig - p(7))) ./ ((1 + exp(-p(5) .* (comSig - p(7)))).^2))];
% 
%                     
%                     
% 
% options.lb = [0.3, 4.8, 0.025, 0.017, 0.1, 40,   4.1];
% options.ub = [0.5, 5.0, 0.100, 0.018, 100, 70,   4.2];
% 
% 
% options.ipopt.limited_memory_update_type = 'bfgs';
% options.ipopt.hessian_approximation = 'limited-memory';
% % options.ipopt.print_level = 1;
% 
% start = [0.4, 4.9, 0.038, 0.0176, 20, 50.53,4.15];
% 
% options.ipopt.max_iter = 100; % max. number of iterations
% 
% % [fitResult, ~] = ipopt (start, funcs, options);
% 
% fitResult = start;
% 
% %%
% close all
% start = fitResult;
% figure
% plot(fitFunction(eneT, sigmasT, fitResult(1), fitResult(2), fitResult(3), fitResult(4), fitResult(5), fitResult(6), fitResult(7)));
% hold on 
% plot(bestT, 'or')
% 
% 
% figure
% plot(fitFunction(eneG, sigmasG, fitResult(1), fitResult(2), fitResult(3), fitResult(4), fitResult(5), fitResult(6), fitResult(7)));
% hold on
% plot(bestG, 'or')
% 
% %[0.450000000000000,4.92065066915907,0.0740426814616639,0.0176945743884910,20.8762496314806,50.5347210423522,4.15591373574955]
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% 
% result = @(E, sigma) fitFunction(E, sigma , 0.4, 4.9, 0.038, 0.0176, 20, 50.53, 4.15);
% 
% countI = 1;
% countJ = 1;
% test = zeros(100,100);
% for i = linspace(50,220,100)
%     countJ = 1;
%     for j = linspace(0,5,100)
%         test(countI,countJ) = result(i,j);
%         countJ = countJ + 1;
%     end
%     countI = countI + 1;
% end
% 
% 
% 
% 
% 
