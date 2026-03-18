function [ct] = helper_createTestCt(cubeDim, resolution, varargin)
% helper_createTestCt  Creates a minimal matRad ct struct for testing.
%
% call
%   [ct] = helper_createTestCt()
%   [ct] = helper_createTestCt(cubeDim)
%   [ct] = helper_createTestCt(cubeDim, resolution)
%
% inputs
%   cubeDim    [nRows nCols nSlices], default [10 12 8]
%              (unequal dimensions help catch x/y coordinate-swap bugs)
%   resolution scalar (isotropic mm), 1x3 vector [x y z mm], or struct
%              with fields .x .y .z; default struct(x=1, y=2, z=3)
%              (unequal resolutions help catch axis-scaling bugs)
%  Options (name-value pairs)
%   'createCoordinateArrays' (logical) if true, adds ct.x ct.y ct
%              ct.z arrays in world mm coordinates; default false
%   'datatype' (string) numeric class for ct.cubeHU; default 'double
%   'HUvalue' (numeric scalar) value to fill ct.cubeHU with; default 0
%
% outputs
%   ct   matRad ct struct with cubeDim and resolution populated

p = inputParser;

p.addOptional('cubeDim', [10 12 8], @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 3, 'positive', 'integer'}));
p.addOptional('resolution', [1 2 3], @(x) (isstruct(x) && all(isfield(x, {'x', 'y', 'z'}))) || ...
              (isnumeric(x) && (isscalar(x) || (isvector(x) && numel(x) == 3))));
p.addParameter('createCoordinateArrays', false, @(x) islogical(x) && isscalar(x));
p.addParameter('datatype', 'double', @(x) ischar(x) || (isstring(x) && isscalar(x)));
p.addParameter('HUvalue', 0, @(x) isnumeric(x) && isscalar(x));

args = varargin;
if nargin >= 2
    args = [{resolution}, args];
end

if nargin >= 1
    args = [{cubeDim}, args];
end

p.parse(args{:});

ct.cubeDim = p.Results.cubeDim;

if isstruct(p.Results.resolution)
    ct.resolution = p.Results.resolution;
elseif isscalar(p.Results.resolution)
    ct.resolution.x = p.Results.resolution;
    ct.resolution.y = p.Results.resolution;
    ct.resolution.z = p.Results.resolution;
else
    ct.resolution.x = p.Results.resolution(1);
    ct.resolution.y = p.Results.resolution(2);
    ct.resolution.z = p.Results.resolution(3);
end

if p.Results.createCoordinateArrays
    ct = matRad_getWorldAxes(ct);
end

ct.cubeHU{1} = p.Results.HUvalue * ones(ct.cubeDim, p.Results.datatype);

end
