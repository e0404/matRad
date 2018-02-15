function matRad_dispToConsole(string,param,typeOfMessage,formatSpec)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad_dispToConsole function to display information on the console 
% depending on the current param.logLevel

% call
%   matRad_dispToConsole(string,param,typeOfMessage)
%
% input
%   string:          message that should be displayed on the console
%   param:           structure containt a subfield logLevel
%   typeOfMessage:   possible options are 'info','warning','error'
%   formatSpec:      optional
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

if exist('param','var')
   if ~isfield(param,'logLevel')
      param.logLevel = 1;
   end
else
   param.logLevel = 1;
end

switch typeOfMessage

   case {'error'}

        fprintf(string);


   case{'warning'}

      if param.logLevel < 4

         fprintf(['[\bWarning:]\b ' string]);

      end

   case{'info'}

       if param.logLevel < 3

          if exist('formatSpec','var')
              fprintf(formatSpec,string);
          else
              fprintf(string);
          end


       end
   otherwise
      
      error('type of message not defined');
      
end
   
   
end