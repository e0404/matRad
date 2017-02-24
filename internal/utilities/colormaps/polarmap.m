function cmap = polarmap(varargin)
%POLARMAP Polarized color map
%	POLARMAP applies a "polarized" blue-white-red colormap to current figure,
%	and adjusts the color axis limits to be centered to zero.
%
%	POLARMAP(M) fixes the number of colors to M (default is 64).
%
%	POLARMAP(MAP) applies linear shading to white to the center of colormap
%	MAP which can be any of existing colormaps (an Mx3 matrix of RGB).
%
%	POLARMAP(MAP,C) uses exponent C to modify the shading contrast. Default 
%	is C = 1 for linear shading. Use C = 2 to strengthen the shading, or 
%	C = 0.5 to attenuate it.
%
%	C=POLARMAP(...) returns an M-by-3 matrix containing the colormap, that 
%	can be used with COLORMAP function like other colormaps.
%
%	Examples:
%		pcolor(peaks), shading interp
%		polarmap, colorbar
%
%	then try the following
%		polarmap(jet,0.5)
%
%	Note the polar shading has no real interest with colormaps that include
%	white color as one of the extremes (like GRAY, BONE, HOT, ...).
%
%	See also JET, HSV, COPPER, SPRING, SUMMER, WINTER, COOL, COLORMAP, RGBPLOT.
%
%	Author: Francois Beauducel, IPGP
%	Created: 2011-10-26
%	Updated: 2012-06-12

%	Copyright (c) 2012, François Beauducel, covered by BSD License.
%	All rights reserved.
%
%	Redistribution and use in source and binary forms, with or without 
%	modification, are permitted provided that the following conditions are 
%	met:
%
%	   * Redistributions of source code must retain the above copyright 
%	     notice, this list of conditions and the following disclaimer.
%	   * Redistributions in binary form must reproduce the above copyright 
%	     notice, this list of conditions and the following disclaimer in 
%	     the documentation and/or other materials provided with the distribution
%	                           
%	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%	POSSIBILITY OF SUCH DAMAGE.

% default parameters
m = 100;	% number of colors
c = 1;	% exponent of shading factor (1 = linear)

if nargin > 0
	if ~isnumeric(varargin{1}) | (size(varargin{1},2) ~= 3 & ~isscalar(varargin{1}))
		error('First argument must be numeric: scalar M or Mx3 color matrix');
	end
	if isscalar(varargin{1})
		m = varargin{1};
	end
end
if nargin > 0 & size(varargin{1},2) == 3
		map = varargin{1};
		m = size(map,1);
else
	map = bluered(m);
end

if nargin > 1 & isscalar(varargin{2})
	c = varargin{2};
end

% linear shading from min/max (colormap value) to center (white)
r = repmat(abs(linspace(1,-1,m)).^c,[3,1])';
map = map.*r + 1 - r;

if nargout > 0
	cmap = map;
else
	colormap(gca,map)
	caxis([-1,1]*max(abs(caxis)))
	% Note: this fixes color axis to manual mode...
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map = bluered(m)

if mod(m,2)
	z = [0,0,0];
	m2 = floor(m/2);
else
	z = zeros([0,3]);
	m2 = m/2;
end
map = [repmat([0,0,1],[m2,1]);z;repmat([1,0,0],[m2,1])];

