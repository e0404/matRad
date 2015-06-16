function [f, g] = matRad_bioObjFunc(w,dij,cst)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad optimization function based on the biologcial effect
% 
% call
%   [f, g] = matRad_bioObjFunc(w,dij,cst)
%
% input
%   w:   weight vector
%   dij: matRad dij struct
%   cst: cst file
%
% output
%   f: objective function value
%   g: gradient vector
%
% References
%   [1] http://iopscience.iop.org/0031-9155/51/12/009
%
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


% calculate biological effect
linTerm  = dij.mAlphaDose*w;
quadTerm = dij.mSqrtBetaDose*w;
e = linTerm + quadTerm.^2;

%Dcut marks the transition from linear quadratic to purely linear shape at
%high doses
d = dij.physicalDose*w;
CutIdx = d>dij.Dcut;

if sum(CutIdx>0)

    linTermHighDose  = (linTerm(CutIdx)./d(CutIdx)).*dij.Dcut;
    quadTermHighDose = (quadTerm(CutIdx)./d(CutIdx)).*dij.Dcut;
    % correct bio effect for voxels having a higher dose than Dcut
    e(CutIdx) = linTermHighDose + quadTermHighDose.^2 +...
        (d(CutIdx)-dij.Dcut).*dij.Smax(CutIdx);
end

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

% Compute optimization function for every VOI.
for  i = 1:size(cst,1)
    
    % Only take OAR or target VOI.
    if ~isempty(cst{i,4}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
        
        % get effect vector in current VOI
        e_i = e(cst{i,4});
        
        % loop over the number of constraints for the current VOI
        for j = 1:size(cst{i,6},1)
            
            % get Penalty
            rho = cst{i,6}(j).parameter(1);
            
            % refernce effect
            if ~isequal(cst{i,6}(j).type, 'mean') && ~isequal(cst{i,6}(j).type, 'EUD')
                e_ref = dij.ax(cst{i,4}).*cst{i,6}(j).parameter(2)+dij.bx(cst{i,4})*cst{i,6}(j).parameter(2)^2;
            end
            
            if isequal(cst{i,6}(j).type, 'square underdosing')
  
                if ~isequal(cst{i,3},'OAR')
                    % underdose : effect minus reference effect
                    underdose = e_i - e_ref;

                    % apply positive operator
                    underdose(underdose>0) = 0;

                    % calculate objective function
                    f = f + (rho/size(cst{i,4},1))*(underdose'*underdose);

                    % calculate delta
                    delta_underdose(cst{i,4}) = delta_underdose(cst{i,4}) + (rho/size(cst{i,4},1))*underdose;
                else
                    disp(['square underdosing constraint for ' cst{i,2} ' will be skipped'])
                end
                
            elseif isequal(cst{i,6}(j).type, 'square overdosing')
                
                % overdose : Dose minus prefered dose
                overdose = e_i - e_ref;
                
                % apply positive operator
                overdose(overdose<0) = 0;
                
                % calculate objective function
                f = f + (rho/size(cst{i,4},1))*(overdose'*overdose);
                
                %calculate delta
                delta_overdose(cst{i,4}) = delta_overdose(cst{i,4}) + (rho/size(cst{i,4},1))*overdose;
                
            elseif isequal(cst{i,6}(j).type, 'square deviation')
                
                if ~isequal(cst{i,3},'OAR')
                    % deviation : Dose minus prefered dose
                    deviation = e_i - e_ref;

                    % claculate objective function
                    f = f + (rho/size(cst{i,4},1))*(deviation'*deviation);

                    % calculate delta
                    delta_deviation(cst{i,4}) = delta_deviation(cst{i,4}) + (rho/size(cst{i,4},1))*deviation;
                else
                    disp(['square deviation constraint for ' cst{i,2} ' will be skipped'])
                end
            
                
            elseif isequal(cst{i,6}(j).type, 'mean')              
                
                if ~isequal(cst{i,3},'TARGET')
                    % calculate objective function
                    f = f + (rho/size(cst{i,4},1))*sum(e_i);

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
                    if sum(e_i.^exponent)>0

                        f = f + rho*nthroot((1/size(cst{i,4},1))*sum(e_i.^exponent),exponent);

                        delta_EUD(cst{i,4}) = delta_EUD(cst{i,4}) + ...
                            rho*nthroot(1/size(cst{i,4},1),exponent) * sum(e_i.^exponent)^((1-exponent)/exponent) * (e_i.^(exponent-1));                    

                        if sum(~isfinite(delta_EUD)) > 0 % check for inf and nan for numerical stability
                            error(['EUD computation for ' cst{i,2} ' failed. Reduce exponent to resolve numerical problems.']);
                        end

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

% gradient calculation
if nargout > 1
 
    delta = 2*(delta_underdose + delta_overdose + delta_deviation) + delta_mean + delta_EUD;        
    vBias= (delta' * dij.mAlphaDose)';
    mPsi = (2*(delta.*quadTerm)'*dij.mSqrtBetaDose)';
    g    =  vBias+mPsi ; 
    
end

