function matRad_dispToConsole(string,typeOfMessage)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_dispToConsole function to display information on the console 
% depending on the current param.logLevel

% call
%   matRad_dispToConsole(string,param,typeOfMessage)
%
% input
%   string:          message that should be displayed on the console
%   typeOfMessage:   possible options are 'info','warning','error'
%
% output
%   
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2017 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global LogLevel 

if ~exist('LogLevel','var')
   LogLevel = 1;
end

switch typeOfMessage

   case {'error'}

        fprintf(string);


   case{'warning'}

      if param.logLevel < 4

         fprintf(string);

      end

   case{'info'}

       if param.logLevel < 3

          fprintf(string);

      end
end
   
   
end