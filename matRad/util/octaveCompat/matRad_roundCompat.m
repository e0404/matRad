function y = matRad_roundCompat(x, n)
% matRad function mimicking Matlab's round(x,n) for compatibility with
% Octave, whose round takes the value to be rounded only and errors with
% "Invalid call to round" as soon as a digit count is passed
%
% call:
%   y = matRad_roundCompat(x)
%   y = matRad_roundCompat(x,n)
%
% input:
%   x   value(s) to round
%   n   number of decimal digits to keep. Positive n rounds to the right of
%       the decimal point, negative n to the left, n = 0 rounds to integers.
%       Omitting n is the same as passing 0
%
% output:
%   y   x rounded to n decimal digits
%
% Note that the Octave branch below scales, rounds and scales back, which is
% not bit-identical to Matlab's round(x,n) - the latter is more accurate for
% values whose scaled magnitude approaches the limit of exact integer
% representation in a double (2^53). Several hot paths in matRad therefore
% keep their own local variant of the scaling idiom instead of calling this
% function, so that they compute the same thing under both environments
% (see matRad_calcBixelWeightAndGradient or matRad_daoVec2ApertureInfo).
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2026 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 || isempty(n)
    y = round(x);
    return
end

% cached, since this is called from loops and the environment cannot change
% within a session
persistent hasDigitArgument
if isempty(hasDigitArgument)
    matRad_cfg = MatRad_Config.instance();
    hasDigitArgument = ~matRad_cfg.isOctave;
end

if hasDigitArgument
    y = round(x, n);
else
    scale = 10.^n;
    y = round(x .* scale) ./ scale;
end

end
