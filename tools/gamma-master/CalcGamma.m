function gamma = CalcGamma(varargin)
% CalcGamma computes 1-D, 2-D, or 3-D global or absolute gamma between two
% datasets (reference and target) given a defined coordinate space. The 
% datasets must have the same number of dimensions, although they can be 
% different sizes. Gamma will be computed for each target dose point by
% shifting the reference image (using linear interpolation) and determining
% the minimum Gamma index across all shifts.
%
% This function optionally uses the Parallel Computing Toolbox GPU interp
% functions to increase computation speed. A try-catch statement is used
% to test for GPU support. In addition, for memory management, the
% meshgrid and data arrays are converted to single precision during
% interpolation. This function calls Event() to log execution status, if 
% available.
%
% For more information on the Gamma evaluation function, see D. A. Low et 
% al., "A technique for the quantitative evaluation of dose distributions", 
% Med Phys. 1998 May; 25(5): 656-61.
%
% The following variables are required for proper execution: 
%   varargin{1}: structure containing the reference data, where the field
%       start is an array containing the coordinates along each dimension
%       of the first voxel, width is an array containing the width of each
%       voxel along each dimension, and data is an n-dimensional array
%   varargin{2}: structure containing the target data, where the field
%       start is an array containing the coordinates along each dimension
%       of the first voxel, width is an array containing the width of each
%       voxel along each dimension, and data is an n-dimensional array
%   varargin{3}: Gamma absolute criterion percentage
%   varargin{4}: Gamma Distance To Agreement (DTA) criterion, in the same
%       units as the reference and target width structure fields  
%   varargin{5:end} (optional): additional parameters preceded by option
%       flags.  The available options are 'local', 'refval', 'restrict',
%       'res', and 'limit'.
%
% The following options can be passed to this argument as name/value pairs:
%   local: boolean, indicates whether to perform a local (1) or global (0) 
%       Gamma computation. If not present, the function will assume a 
%       global Gamma computation.
%   refval: reference value for the global absolute criterion. Is used with 
%       the percentage from varargin{3} to compute absolute value. If not 
%       present, the maximum value in the reference data is used.
%   restrict: restricted search flag. If 1, only the gamma values along the 
%       X/Y/Z axes are computed during 3D comptation. If 0 or not provided, 
%       the entire rectangular search space is computed.
%   res: DTA resolution factor.  The number of distance steps equal the
%       resolution factor multiplied by the limit.  If not provided, the
%       factor is 100 for 1D/2D and 20 for 3D calculations.
%   limit: The DTA limit.  This number determines how far the function will 
%       search in the distance axes when computing Gamma.  This also 
%       therefore specifies the maximum "real" Gamma Index value. 
%   cpu: boolean, set to 1 to force CPU computation.
%
% The following variables are returned upon succesful completion:
%   gamma: array of the same dimensions as varargin{2}.data containing the
%       computed gamma values
%
% Below is an example of how the function is used:
%
%   reference.start = [-10 -10]; % mm
%   reference.width = [0.1 0.1]; % mm
%   reference.data = rand(200);
%
%   target.start = [-10 -10]; % mm
%   target.width = [0.1 0.1]; % mm
%   target.data = rand(200);
%
%   percent = 3;
%   dta = 0.5; % mm
%   local = 0; % Perform global gamma
%   
%   gamma = CalcGamma(reference, target, percent, dta, 'local', local);
%
% Author: Mark Geurts, mark.w.geurts@gmail.com
% Copyright (C) 2014 University of Wisconsin Board of Regents
%
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the  
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along 
% with this program. If not, see http://www.gnu.org/licenses/.

%% Validate Inputs
% Verify at least four input arguments are provided
if nargin < 4
    
    % If not, throw an error and stop execution
    if exist('Event', 'file') == 2
        Event('Too few input argumenst passed to CalcGamma', 'ERROR');
    else
        error('Too few input argumenst passed to CalcGamma');
    end
end

% Check if the reference structure contains width, start, and data fields,
% and if the size of the width and start vectors are equal
if ~isfield(varargin{1}, 'width') || ~isfield(varargin{1}, 'start') || ...
        ~isfield(varargin{1}, 'data') || ~isequal(size(varargin{1}.width), ...
        size(varargin{1}.start))
    
    % If not, throw an error and stop execution
    if exist('Event', 'file') == 2
        Event(['Incorrect reference data format.  Must contain width, ', ...
            'start, and data fields and be of equal dimensions'], 'ERROR');
    else
        error(['Incorrect reference data format.  Must contain width, ', ...
            'start, and data fields and be of equal dimensions']);
    end
    
% Check if the target structure contains width, start, and data fields,
% and if the size of the width and start vectors are equal
elseif ~isfield(varargin{2}, 'width') || ~isfield(varargin{2}, 'start') || ...
        ~isfield(varargin{2}, 'data') || ~isequal(size(varargin{2}.width), ...
        size(varargin{2}.start))
    
    % If not, throw an error and stop execution
    if exist('Event', 'file') == 2
        Event(['Incorrect target data format.  Must contain width, ', ...
            'start, and data fields and be of equal dimensions'], 'ERROR');
    else
        error(['Incorrect target data format.  Must contain width, ', ...
            'start, and data fields and be of equal dimensions']);
    end
    
% Check if the reference and target data arrays are the same number of
% dimensions.  Calculating the gamma from a lower dimensional dataset to a
% higher dimensional reference is currently not supported
elseif ~isequal(size(size(varargin{1}.data)), size(size(varargin{2}.data)))
    
    % If not, throw an error and stop execution
    if exist('Event', 'file') == 2
        Event(['The fixed and target data arrays must be of the same', ...
            ' dimensions'], 'ERROR');
    else
        error(['The fixed and target data arrays must be of the same', ...
            ' dimensions']);
    end
end

% Log validation completed
if exist('Event', 'file') == 2
    Event('Data validation completed, beginning Gamma computation');
    tic;
end

%% Declare default options
% local indicates whether to perform a local (1) or global (0) Gamma 
% computation. If not present, the function will assume a global Gamma 
% computation.
local = 0;

% refval is the reference value for the global absolute criterion. Is used 
% with the percentage from varargin{3} to compute absolute value. If not 
% present, the maximum value in the reference data is used.
refval = max(max(max(varargin{1}.data)));

% restricted search flag. If 1, only the gamma values along the X/Y/Z axes 
% are computed during 3D comptation. If 0 or not provided, the entire 
% rectangular search space is computed.
restrict = 0;

% The resolution parameter determines the number of steps (relative to 
% the distance to agreement) that each reference voxel will be
% interpolated to and gamma calculated.  A value of 5 with a DTA of 3
% mm means that gamma will be calculated at intervals of 3/5 = 0.6 mm.
% Different resolutions can be set for different dimensions of data. 
if size(varargin{2}.width,2) == 1

    % Set 1-D resolution
    res = 100;
    
elseif size(varargin{2}.width,2) == 2

    % Set 2-D resolution
    res = 50;
    
elseif size(varargin{2}.width,2) == 3

    % Set 3-D resolution
    res = 20;
    
end

% The search limit parameter determines how far the function will search in
% the distance axes when computing Gamma.  This also therefore specifies
% the maximum believable Gamma Index value.  Typically this value is 2.
% This should always be set to an integer.
limit = 2;


%% Search for provided options
% Load data structure from varargin
for i = 5:nargin
    
    % If the local option is set
    if strcmp(varargin{i}, 'local')
        local = varargin{i+1}; 
        
    % If the refval option is set
    elseif strcmp(varargin{i}, 'refval')
        refval = varargin{i+1}; 
        
    % If the restrict option is set
    elseif strcmp(varargin{i}, 'restrict')
        restrict = varargin{i+1}; 
        
    % If the res option is set
    elseif strcmp(varargin{i}, 'res')
        res = varargin{i+1}; 
        
    % If the limit option is set
    elseif strcmp(varargin{i}, 'limit')
        limit = varargin{i+1}; 
    
    % If the cpu option is set
    elseif strcmp(varargin{i}, 'cpu')
        cpu = varargin{i+1}; 
    end
end

%% Log options
% If Event reporting is enabled
if exist('Event', 'file') == 2
    
    % Log local
    if local == 1
        Event('Gamma calculation set to local');  
    else
        Event('Gamma calculation assumed to global');
    end
    
    % Log refval
    Event(sprintf('Reference value set to %g', refval));
    
    % Log restrict
    if restrict == 1
        Event('Restricted search enabled');  
    else
        Event('Restricted search disabled');
    end

    % Log res
    Event(sprintf('Resolution set to %g', res));
    
    % Log limit
    Event(sprintf('DTA limit set to %g', limit));
end

%% Compute mesh grids
% If the reference dataset is 1-D
if size(varargin{1}.width,2) == 1

    % Log event
    if exist('Event', 'file') == 2
        Event('Reference dataset is 1-D');
    end
    
    % Check if the data is in rows or columns (this is only needed for 1-D)
    if size(varargin{1}.data,1) > size(varargin{1}.data,2)
        
        % If in rows, transpose
        varargin{1}.data = varargin{1}.data';
    end
    
    % Compute the reference X coordinates using the start and width values
    refX = single(varargin{1}.start(1):varargin{1}.width(1):varargin{1}.start(1) ...
        + varargin{1}.width(1) * (size(varargin{1}.data,2) - 1));
    
% Otherwise, if the reference dataset is 2-D
elseif size(varargin{1}.width,2) == 2

    % Log event
    if exist('Event', 'file') == 2
        Event('Reference dataset is 2-D');
    end
    
    % Compute X and Y meshgrids for the reference dataset positions using 
    % the start and width values
    [refX, refY] = meshgrid(single(varargin{1}.start(2): ...
        varargin{1}.width(2):varargin{1}.start(2) + varargin{1}.width(2) * ...
        (size(varargin{1}.data,2) - 1)), single(varargin{1}.start(1): ...
        varargin{1}.width(1):varargin{1}.start(1) + varargin{1}.width(1)...
        * (size(varargin{1}.data,1) - 1)));
    
% Otherwise, if the reference dataset is 3-D
elseif size(varargin{1}.width,2) == 3

    % Log event
    if exist('Event', 'file') == 2
        Event('Reference dataset is 3-D');
    end
    
    % Compute X, Y, and Z meshgrids for the reference dataset positions
    % using the start and width values, permuting X/Y
    [refX, refY, refZ] = meshgrid(single(varargin{1}.start(2): ...
        varargin{1}.width(2):varargin{1}.start(2) + varargin{1}.width(2) * ...
        (size(varargin{1}.data,2) - 1)), single(varargin{1}.start(1): ...
        varargin{1}.width(1):varargin{1}.start(1) + varargin{1}.width(1)...
        * (size(varargin{1}.data,1) - 1)), single(varargin{1}.start(3):...
        varargin{1}.width(3):varargin{1}.start(3) + varargin{1}.width(3)...
        * (size(varargin{1}.data,3) - 1)));

% Otherwise, if the reference data is of higher dimension
else

    % Throw an error and stop execution
    if exist('Event', 'file') == 2
        Event('The fixed data structure contains too many dimensions', ...
            'ERROR');
    else
        error('The fixed data structure contains too many dimensions');
    end
end

% If the target dataset is 1-D
if size(varargin{2}.width,2) == 1

    % Log event
    if exist('Event', 'file') == 2
        Event('Target dataset is 1-D');
    end
    
    % Check if the data is in rows or columns (this is only needed for 1-D)
    if size(varargin{2}.data,1) > size(varargin{2}.data,2)
        
        % If in rows, transpose
        varargin{2}.data = varargin{2}.data';
    end
    
    % Compute the target X coordinates using the start and width values
    tarX = single(varargin{2}.start(1):varargin{2}.width(1):varargin{2}.start(1) ...
        + varargin{2}.width(1) * (size(varargin{2}.data,2) - 1));
    
% Otherwise, if the target dataset is 2-D
elseif size(varargin{2}.width,2) == 2

    % Log event
    if exist('Event', 'file') == 2
        Event('Target dataset is 2-D');
    end
    
    % Compute X and Y meshgrids for the target dataset positions using the
    % start and width values
    [tarX, tarY] = meshgrid(single(varargin{2}.start(2):...
        varargin{2}.width(2):varargin{2}.start(2) + varargin{2}.width(2) * ...
        (size(varargin{2}.data,2) - 1)), single(varargin{2}.start(1): ...
        varargin{2}.width(1):varargin{2}.start(1) + varargin{2}.width(1) ...
        * (size(varargin{2}.data,1) - 1)));
    
% Otherwise, if the target dataset is 3-D
elseif size(varargin{2}.width,2) == 3

    % Log event
    if exist('Event', 'file') == 2
        Event('Target dataset is 3-D');
    end
    
    % Compute X, Y, and Z meshgrids for the target dataset positions using
    % the start and width values, permuting X/Y
    [tarX, tarY, tarZ] = meshgrid(single(varargin{2}.start(2):...
        varargin{2}.width(2):varargin{2}.start(2) + varargin{2}.width(2) * ...
        (size(varargin{2}.data,2) - 1)), single(varargin{2}.start(1): ...
        varargin{2}.width(1):varargin{2}.start(1) + varargin{2}.width(1) ...
        * (size(varargin{2}.data,1) - 1)), single(varargin{2}.start(3):...
        varargin{2}.width(3):varargin{2}.start(3) + varargin{2}.width(3) ...
        * (size(varargin{2}.data,3) - 1)));
    
end

%% Initialize variables
% Generate an initial gamma volume with values of 2^2 (this is the maximum
% reliable value of gamma, see description of limit above).  Note that
% gamma-squared is stored during computation; sqrt is computed at the end.
gamma = ones(size(varargin{2}.data)) * (limit ^ 2);

% Log number of gamma calculations (for status updates on 3D calcs)
if restrict == 1

    % Compute number of restricted search calcs
    num = res * (limit * 2) * size(varargin{2}.width, 2);
else

    % Compute total number of calcs
    num = res * (limit * 2) ^ size(varargin{2}.width, 2);
end

% num is the number of iterations, num * numel the total number of
% interpolations being performed
if exist('Event', 'file') == 2
    Event(sprintf('Number of gamma calculations = %g', num * numel(gamma)));
end

% Initialize counter (for progress indicator)
n = 0;

%% Begin computation
% Start try-catch block to safely test for CUDA functionality
try
    % If the cpu option is set, throw an error to force CPU computation
    if exist('cpu', 'var') == 1 && cpu == 1
        error('Reverting to CPU computation');
    end
    
    % Clear and initialize GPU memory.  If CUDA is not enabled, or if the
    % Parallel Computing Toolbox is not installed, this will error, and the
    % function will automatically rever to CPU computation via the catch
    % statement
    gpuDevice(1);
    
    % Start a for loop to interpolate the dose array along the x-direction.  
    % Note to support parfor loops indices must be integers, so x varies 
    % from -2 to +2 multiplied by the number of interpolation steps.  
    % Effectively, this evaluates gamma from -2 * DTA to +2 * DTA.
    for x = -limit*res:limit*res
        
        % i is the x axis step value
        i = x/res * varargin{4};
        
        % Initialize j and k as zero (they will be updated if the data is
        % of higher dimension)
        j = 0;
        k = 0;
        
        % If the data contains a second dimension
        if size(varargin{1}.width,2) > 1
   
            % Start a for loop to interpolate the dose array along the
            % y-direction.  Note to support parfor loops indices must be
            % integers, so y varies from -2 to +2 multiplied by the number
            % of interpolation steps.  Effectively, this evaluates gamma
            % from -2 * DTA to +2 * DTA.
            for y = -limit*res:limit*res
                
                % j is the y axis step value
                j = y/res * varargin{4};
                
                % Initialize k as zero (it will be updated if the data is
                % of higher dimension)
                k = 0;
                
                % If the data contains a third dimension
                if size(varargin{1}.width, 2) > 2
                    
                    % Start a for loop to interpolate the dose array along 
                    % the z-direction.  Note to support parfor loops 
                    % indices must be integers, so z varies from -2 to +2 
                    % multiplied by the number of interpolation steps.
                    % Effectively, this evaluates gamma from -2 * DTA to 
                    % +2 * DTA.
                    for z = -limit*res:limit*res
                        
                        % k is the z axis step value
                        k = z/res * varargin{4};

                        % Check restricted search flag
                        if restrict == 0 || sum(abs([x y z]) > 0) == 1
                            
                            % Run GPU interp3 function to compute the reference
                            % values at the specified target coordinate points
                            interp = gather(interp3(gpuArray(refX), gpuArray(refY), ...
                                gpuArray(refZ), gpuArray(single(varargin{1}.data)), ...
                                gpuArray(tarX + i), gpuArray(tarY + j), ...
                                gpuArray(tarZ + k), 'linear', 0));

                            % Update the gamma array by returning the minimum
                            % of the existing value or the new value
                            gamma = min(gamma, GammaEquation(interp, ...
                                varargin{2}.data, i, j, k, varargin{3}, varargin{4}, ...
                                refval, local));
                            
                            % Update counter 
                            n = n + 1;
                            
                            % If counter is at an even %, display progress
                            if mod((n-1)/num, 0.01) > 0.005 && ...
                                    mod(n/num, 0.01) < 0.005
                                fprintf('%0.1f%%\n', n/num*100);
                            end
                        end
                    end
                    
                % Otherwise, the data is 2-D
                else
                
                    % Run GPU interp2 function to compute the reference
                    % values at the specified target coordinate points
                    interp = gather(interp2(gpuArray(refX), gpuArray(refY), ...
                        gpuArray(single(varargin{1}.data)), gpuArray(tarX + i), ...
                        gpuArray(tarY + j), 'linear', 0));
                    
                    % Update the gamma array by returning the minimum
                    % of the existing value or the new value
                    gamma = min(gamma, GammaEquation(interp, varargin{2}.data, ...
                        i, j, k, varargin{3}, varargin{4}, ...
                        refval, local));
                end
            end
            
        % Otherwise, the data is 1-D
        else
        
            % Run GPU interp function to compute the reference values at 
            % the specified target coordinate points
            interp = gather(interp1(gpuArray(refX), ...
                gpuArray(single(varargin{1}.data)), gpuArray(tarX + i), ...
                'linear', 0));
            
            % Update the gamma array by returning the minimum of the 
            % existing value or the new value
            gamma = min(gamma, GammaEquation(interp, varargin{2}.data, ...
                i, j, k, varargin{3}, varargin{4}, ...
                refval, local));
        end
    end
    
% If GPU fails, revert to CPU computation
catch

    % Log GPU failure (if cpu flag is not set)
    if exist('Event', 'file') == 2 && exist('cpu', 'var') ~= 1
        Event('GPU failed, reverting to CPU method', 'WARN'); 
    end
    
    % Start a for loop to interpolate the dose array along the x-direction.  
    % Note to support parfor loops indices must be integers, so x varies 
    % from -2 to +2 multiplied by the number of interpolation steps.  
    % Effectively, this evaluates gamma from -2 * DTA to +2 * DTA.
    for x = -limit*res:limit*res
    
        % i is the x axis step value
        i = x/res * varargin{4};
        
        % Initialize j and k as zero (they will be updated if the data is
        % of higher dimension)
        j = 0;
        k = 0;
        
        % If the data contains a second dimension
        if size(varargin{1}.width, 2) > 1
            
            % Start a for loop to interpolate the dose array along the
            % y-direction.  Note to support parfor loops indices must be
            % integers, so y varies from -2 to +2 multiplied by the number
            % of interpolation steps.  Effectively, this evaluates gamma
            % from -2 * DTA to +2 * DTA.
            for y = -limit*res:limit*res
                
                % j is the y axis step value
                j = y/res * varargin{4};
                
                % Initialize k as zero (it will be updated if the data is
                % of higher dimension)
                k = 0;
                
                % If the data contains a third dimension
                if size(varargin{1}.width, 2) > 2
                    
                    % Start a for loop to interpolate the dose array along 
                    % the z-direction.  Note to support parfor loops 
                    % indices must be integers, so z varies from -2 to +2 
                    % multiplied by the number of interpolation steps.
                    % Effectively, this evaluates gamma from -2 * DTA to 
                    % +2 * DTA.
                    for z = -limit*res:limit*res
                        
                        % k is the z axis step value
                        k = z/res * varargin{4};

                        % Check restricted search flag
                        if restrict == 0 || sum(abs([x y z]) > 0) == 1
                            
                            % Run CPU interp3 function to compute the reference
                            % values at the specified target coordinate points
                            interp = interp3(refX, refY, refZ, ...
                                single(varargin{1}.data), tarX + i, ...
                                tarY + j, tarZ + k, '*linear', 0);

                            % Update the gamma array by returning the minimum
                            % of the existing value or the new value
                            gamma = min(gamma, GammaEquation(interp, ...
                                varargin{2}.data, i, j, k, varargin{3}, ...
                                varargin{4}, refval, local));
                            
                            % Update counter 
                            n = n + 1;
                            
                            % If counter is at an even %, display progress
                            if mod((n-1)/num, 0.01) > 0.005 && ...
                                    mod(n/num, 0.01) < 0.005
                                fprintf('%0.1f%%\n', n/num*100);
                            end
                        end
                    end
                    
                % Otherwise, the data is 2-D
                else
                
                    % Run CPU interp2 function to compute the reference
                    % values at the specified target coordinate points
                    interp = interp2(refX, refY, single(varargin{1}.data), ...
                        tarX + i, tarY + j, '*linear', 0);
                    
                    % Update the gamma array by returning the minimum
                    % of the existing value or the new value
                    gamma = min(gamma, GammaEquation(interp, ...
                        varargin{2}.data, i, j, k, varargin{3}, ...
                        varargin{4}, refval, local));
                end
            end
            
        % Otherwise, the data is 1-D
        else
        
            % Run CPU interp function to compute the reference values at 
            % the specified target coordinate points
            interp = interp1(refX, single(varargin{1}.data), tarX + i, ...
                '*linear', 0);
            
            % Update the gamma array by returning the minimum of the 
            % existing value or the new value
            gamma = min(gamma, GammaEquation(interp, varargin{2}.data, ...
                i, j, k, varargin{3}, varargin{4}, ...
                refval, local));
        end
    end
end
    
%% Finish up
% Take square root of result
gamma = sqrt(gamma);

% Log completion
if exist('Event', 'file') == 2
    Event(sprintf(['Gamma calculation completed successfully in ', ...
        '%0.3f seconds'], toc));
end

% Clear temporary variables
clear refX refY refZ tarX tarY tarZ interp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gamma = GammaEquation(ref, tar, i, j, k, perc, dta, refval, local)
% GammaEquation is the programmatic form of the Gamma definition as given
% by Low et al in matrix form.  This function computes both local and
% global Gamma, and is a subfunction for CalcGamma.  Note that this
% equation returns Gamma-squared instead of Gamma.
%
% The following inputs are used for computation and are required:
%   ref: the reference 3D array.  Must be the same size as tar
%   tar: the target 3D array.  Must be the same size as ref
%   i: magnitude of x position offset of tar to ref, relative to dta
%   j: magnitude of y position offset of tar to ref, relative to dta
%   k: magnitude of z position offset of tar to ref, relative to dta
%   perc: the percent Gamma criterion, given in % (i.e. 3 for 3%)
%   dta: the distance to agreement Gamma criterion, unitless but relative
%       to i, j, and k
%   refval: if global, the reference value to base the % criterion from 
%   local: boolean, indicates whether to perform a local (1) or global (0)
%       Gamma computation
%
% The following variables are returned:
%   gamma: a 3D array of the same dimensions as ref and interp of the
%       computed gamma-squared value for each voxel based on interp and 
%       i, j, k

% If local is set to 1, perform a local Gamma computation
if local == 1
    
    % Gamma is defined as the sqrt((abs difference/relative tolerance)^2 +
    % sum((voxel offset/dta)^2)).  The sqrt is removed for computational
    % efficiency.
    gamma = ((tar - ref) ./ (ref * perc / 100)).^2 + ...
        (i/dta)^2 + (j/dta)^2 + (k/dta)^2;
else
    
    % Gamma is defined as the sqrt((abs difference/absolute  tolerance)^2 +
    % sum((voxel offset/dta)^2)). The sqrt is removed for computational
    % efficiency.
    gamma = ((tar - ref) ./ (refval * perc / 100)).^2 + ...
        (i/dta)^2 + (j/dta)^2 + (k/dta)^2;
end
