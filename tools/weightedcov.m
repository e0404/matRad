function C = weightedcov(Y, w)
%   Weighted Covariance Matrix
%
%   WEIGHTEDCOV returns a symmetric matrix C of weighted covariances
%   calculated from an input T-by-N matrix Y whose rows are
%   observations and whose columns are variables and an input T-by-1 vector
%   w of weights for the observations. This function may be a valid
%   alternative to COV if observations are not all equally relevant
%   and need to be weighted according to some theoretical hypothesis or
%   knowledge.
%
%   C = WEIGHTEDCOV(Y, w) returns a positive semidefinite matrix C, i.e. all its
%   eigenvalues are non-negative.
%
%   If w = ones(size(Y, 1), 1), no difference exists between
%   WEIGHTEDCOV(Y, w) and COV(Y, 1).
%
%   REFERENCE: mathematical formulas in matrix notation are available in
%   F. Pozzi, T. Di Matteo, T. Aste,
%   "Exponential smoothing weighted correlations",
%   The European Physical Journal B, Volume 85, Issue 6, 2012.
%   DOI:10.1140/epjb/e2012-20697-x. 
%
% % ======================================================================
% % EXAMPLE
% % ======================================================================
%
% % GENERATE CORRELATED STOCHASTIC PROCESSES
%   T = 100;                                                                      % number of observations
%   N = 500;                                                                      % number of variables
%   Y = randn(T, N);                                                              % shocks from standardized normal distribution
%   Y = cumsum(Y);                                                                % correlated stochastic processes
%
% % CHOOSE EXPONENTIAL WEIGHTS
%   alpha = 2 / T;
%   w0 = 1 / sum(exp(((1:T) - T) * alpha));
%   w = w0 * exp(((1:T) - T) * alpha);                                            % weights: exponential decay
%
% % COMPUTE WEIGHTED COVARIANCE MATRIX
%   c = weightedcov(Y, w);                                                        % Weighted Covariance Matrix
%
% % ======================================================================
%
%   See also CORRCOEF, COV, STD, MEAN.
%   Check also WEIGHTEDCORRS (FE 20846) and KENDALLTAU (FE 27361)
%
% % ======================================================================
%
%-*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-*%
%                                                                                               %
%            Author: Liber Eleutherios                                             %
%            E-Mail: libereleutherios@gmail.com                             %
%            Date: 15 June 2012                                                    %
%                                                                                               %
%-*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-*%
%
% % ======================================================================
%

% Check input
ctrl = isvector(w) & isreal(w) & ~any(isnan(w)) & ~any(isinf(w)) & all(w > 0);
if ctrl
  w = w(:) / sum(w);                                                              % w is column vector
else
  error('Check w: it needs be a vector of real positive numbers with no infinite or nan values!')
end
ctrl = isreal(Y) & ~any(isnan(Y)) & ~any(isinf(Y)) & (size(size(Y), 2) == 2);
if ~ctrl
  error('Check Y: it needs be a 2D matrix of real numbers with no infinite or nan values!')
end
ctrl = length(w) == size(Y, 1);
if ~ctrl
  error('size(Y, 1) has to be equal to length(w)!')
end

[T, N] = size(Y);                                                                 % T: number of observations; N: number of variables
C = Y - repmat(w' * Y, T, 1);                                                     % Remove mean (which is, also, weighted)
C = C' * (C .* repmat(w, 1, N));                                                  % Weighted Covariance Matrix
C = 0.5 * (C + C');                                                               % Must be exactly symmetric
