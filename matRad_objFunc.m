function [f, g] = matRad_objFunc(w,dij,cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad objective function for inverse planning supporting mean dose
% objectives, EUD objectives, squared overdosage, squared underdosage, and
% squared deviation
% 
% call
%   [f, g] = matRad_objFunc(w,dij,cst)
%
% input
%   w:   bixel weight vector
%   dij: dose influence matrix
%   cst: matRad cst struct
%
% output
%   f: objective function value
%   g: gradient
%
% References
%   [1] http://www.sciencedirect.com/science/article/pii/S0958394701000577
%   [2] http://www.sciencedirect.com/science/article/pii/S0360301601025858
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015, Mark Bangert, on behalf of the matRad development team
%
% m.bangert@dkfz.de
%
% This file is part of matRad.
%
% matrad is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free 
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% matRad is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License in the
% file license.txt along with matRad. If not, see
% <http://www.gnu.org/licenses/>.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate dose
d = dij.physicalDose*w;

% Numbers of voxels
numVoxels = size(dij.physicalDose,1);

% Initializes f
f = 0;

% Initializes delta
delta_underdose = zeros(numVoxels,1);
delta_overdose  = zeros(numVoxels,1);
delta_deviation = zeros(numVoxels,1);
delta_mean      = zeros(numVoxels,1);
delta_EUD       = zeros(numVoxels,1);

% compute objective function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
    
        % get dose vector in current VOI
        d_i = d(cst{i,4});
                
        % loop over the number of constraints for the current VOI
        for j = 1:size(cst{i,6},1)
            
            % get Penalty
            rho = cst{i,6}(j).parameter(1);
            
            if isequal(cst{i,6}(j).type, 'square underdosing') 
               
                if ~isequal(cst{i,3},'OAR')
                    % underdose : Dose minus prefered dose
                    underdose = d_i - cst{i,6}(j).parameter(2);

                    % apply positive operator
                    underdose(underdose>0) = 0;

                    % calculate objective function
                    f = f + (rho/size(cst{i,4},1))*(underdose'*underdose);

                    % calculate delta
                    delta_underdose(cst{i,4}) = delta_underdose(cst{i,4}) +...
                        (rho/size(cst{i,4},1))*underdose;
                else
                    disp(['square underdosing constraint for ' cst{i,2} ' will be skipped'])
                end
                
            elseif isequal(cst{i,6}(j).type, 'square overdosing')
                
                    % overdose : Dose minus prefered dose
                    overdose = d_i - cst{i,6}(j).parameter(2);

                    % apply positive operator
                    overdose(overdose<0) = 0;

                    % calculate objective function
                    f = f + (rho/size(cst{i,4},1))*(overdose'*overdose);

                    %calculate delta
                    delta_overdose(cst{i,4}) = delta_overdose(cst{i,4}) + ...
                        (rho/size(cst{i,4},1))*overdose;
                
           elseif isequal(cst{i,6}(j).type, 'square deviation')
               
               if ~isequal(cst{i,3},'OAR')
                    % deviation : Dose minus prefered dose
                    deviation = d_i - cst{i,6}(j).parameter(2);

                    % claculate objective function
                    f = f + (rho/size(cst{i,4},1))*(deviation'*deviation);

                    % calculate delta
                    delta_deviation(cst{i,4}) = delta_deviation(cst{i,4}) +...
                        (rho/size(cst{i,4},1))*deviation;
                else
                    disp(['square deviation constraint for ' cst{i,2} ' will be skipped'])
                end
                
            elseif isequal(cst{i,6}(j).type, 'mean')              
                
                if ~isequal(cst{i,3},'TARGET')
                    % calculate objective function
                    f = f + (rho/size(cst{i,4},1))*sum(d_i);

                    % calculate delta
                    delta_mean(cst{i,4}) = delta_mean(cst{i,4}) + ...
                        (rho/size(cst{i,4},1))*ones(size(cst{i,4},1),1);
                else
                    disp(['mean constraint for ' cst{i,2} ' will be skipped'])
                end
                
            elseif isequal(cst{i,6}(j).type, 'EUD') 
               
               if ~isequal(cst{i,3},'TARGET')
                    % get exponent for EUD
                    exponent = cst{i,6}(j).parameter(2);

                    % calculate objective function and delta
                    if sum(d_i.^exponent)>0

                        f = f + rho*nthroot((1/size(cst{i,4},1))*sum(d_i.^exponent),exponent);

                        delta_EUD(cst{i,4}) = delta_EUD(cst{i,4}) + ...
                            rho*nthroot(1/size(cst{i,4},1),exponent) * sum(d_i.^exponent)^((1-exponent)/exponent) * (d_i.^(exponent-1));

                    end

                    if sum(~isfinite(delta_EUD)) > 0 % check for inf and nan for numerical stability
                        error(['EUD computation for ' cst{i,2} ' failed. Reduce exponent to resolve numerical problems.']);
                    end
               else
                    disp(['EUD constraint for ' cst{i,2} ' will be skipped'])
                end
               
            else
                
                error('undefined objective in cst struct');
                
            end
      
        end
        
    end
end
   
if nargout > 1
    % Calculate gradient
    g = (( 2*(delta_underdose + delta_overdose + delta_deviation) + delta_mean + delta_EUD )' * dij.physicalDose)';
end

end