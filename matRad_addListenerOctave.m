% Copyright (C) 2020 nikla
% 
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see
% <https://www.gnu.org/licenses/>.

% -*- texinfo -*- 
% @deftypefn {} {@var{retval} =} matRad_addListenerOctave (@var{input1}, @var{input2})
%
% @seealso{}
% @end deftypefn

% Author: nikla <nikla@LAPTOP-NIKLAS>
% Created: 2020-06-04

function retval = matRad_addListenerOctave (hSource,eventName,callback)
    global eventMap;
    
    persistent warningprinted;
    if isempty(warningprinted)
        hSource.showWarning('You are using an experimental implementation of ''addlistener'' concept for Octave');
        warningprinted=2;
    end
    
    newListener = struct('src',hSource,'event',eventName,'callback',callback);
    if isempty(eventMap)
        eventMap = newListener;
    else
        eventMap(end+1) = newListener;
    end
    
    
end
