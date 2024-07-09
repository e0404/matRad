function hlut = matRad_readHLUT(filename)
% matRad function to read HLUT from filename
%
% call
%   hlut = matRad_readHLUT(filename)
%
% input
%   filename: hlut filename
%
% output
%   hlut:     lookup table
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Check the extension
  [~,~,ext] = fileparts(filename);
  
  try
  switch ext
      %% Iterate through supported extensions
      case '.hlut'
          hlutFile = fopen(filename,'r');
          hlut = cell2mat(textscan(hlutFile,'%f %f','CollectOutput',1,'commentStyle','#'));
          fclose(hlutFile);
          
      otherwise 
            error(['HLUT extension ''' ext ''' not supported!']);
  end
  catch ME
      hlut = [];
      error(getReport(ME));      
  end

end

