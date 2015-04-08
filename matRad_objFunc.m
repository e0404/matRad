function [f, g] = matRad_objFunc(w,dij,cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad objective function for inverse planning
% 
% call
%   [f, g, d] = matRad_IMRTObjFunc(w,dij,cst)
%
% input
%   w:   bixel weight vector
%   dij: dose influence matrix
%   cst: matRad cst struct
%
% output
%   f: objective function value
%   g: gradient
%   d: dose distribution
%
% References
%   -
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
    if isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET')
        
        
        % get current indices of VOI
        idx = cst{i,4};
                
        % loop over the number of constraints for the current VOI
        for j = 1:size(cst{i,6},2)
            
            % get current priority
            prior = cst{i,6}(j).priority;
            
            for k=1:size(cst,1)
             for l=1:size(cst{k,6},2)
                if cst{k,6}(l).priority<prior && ~(k==i&&l==j)
                    % remove indices from VOI with higher priority from
                    % current VOI
                    idx=setdiff(idx,cst{k,4});
                end
             end
            end
            
            % get dose vector in current VOI
            d_i = d(idx);
            % get Penalty
            rho = cst{i,6}(j).parameter(1);
            
            if isequal(cst{i,6}(j).type, 'square underdosing') && ~isempty(d_i)
  
                % underdose : Dose minus prefered dose
                underdose = d_i - cst{i,6}(j).parameter(2);
                
                % apply positive operator
                underdose(underdose>0) = 0;
                
                % calculate objective function
                f = f + (rho/size(idx,1))*(underdose'*underdose);
                
                % calculate delta
                delta_underdose(idx) = delta_underdose(idx) +...
                    (rho/size(idx,1))*underdose;
                
            elseif isequal(cst{i,6}(j).type, 'square overdosing') && ~isempty(d_i)
                
                % overdose : Dose minus prefered dose
                overdose = d_i - cst{i,6}(j).parameter(2);
                
                % apply positive operator
                overdose(overdose<0) = 0;
                
                % calculate objective function
                f = f + (rho/size(idx,1))*(overdose'*overdose);
                
                %calculate delta
                delta_overdose(idx) = delta_overdose(cst{i,4}) + ...
                    (rho/size(idx,1))*overdose;
                
            elseif isequal(cst{i,6}(j).type, 'square deviation') && ~isempty(d_i)
                
                % deviation : Dose minus prefered dose
                deviation = d_i - cst{i,6}(j).parameter(2);
                
                % claculate objective function
                f = f + (rho/size(idx,1))*(deviation'*deviation);
                
                % calculate delta
                delta_deviation(idx) = delta_deviation(idx) +...
                    (rho/size(idx,1))*deviation;
                
            elseif isequal(cst{i,6}(j).type, 'mean') && ~isempty(d_i)              
                
                % calculate objective function
                f = f + (rho/size(cst{i,4},1))*sum(d_i);
                
                % calculate delta
                delta_mean(idx) = delta_mean(idx) + ...
                    (rho/size(idx,1))*ones(size(idx,1),1);
                
            elseif isequal(cst{i,6}(j).type, 'EUD') && ~isempty(d_i)
                
                % get exponent for EUD
                exponent = cst{i,6}(j).exponent;
                
                % calculate objective function and delta
                if sum(d_i.^exponent)>0
                    
                    f = f + rho*nthroot((1/size(idx,1))*sum(d_i.^exponent),exponent);
                    
                    delta_EUD(idx) = delta_EUD(idx) + ...
                        rho*nthroot(1/size(idx,1),exponent) * sum(d_i.^exponent)^((1-exponent)/exponent) * (d_i.^(exponent-1));
                    
                end
                
            elseif ~isempty(d_i)
                
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