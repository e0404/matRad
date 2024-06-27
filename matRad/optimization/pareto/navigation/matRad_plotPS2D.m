function matRad_plotPS2D(retStruct,idx)
    % matRad_plotParetoFront3D implements a function to visualize a 3d Pareto
    % surface
    % 
    %
    % input
    %   retStruct:  Structure returned from Pareto optimization
    %   idx(optional):
    % output
    %   resultGUI:  struct containing optimized fluence vector, dose, and (for
    %               biological optimization) RBE-weighted dose etc.
    %   optimizer:  Used Optimizer Object

    %
    % References
    %   
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2024 the matRad development team. 
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ps = retStruct.finds;
    if size(ps,2) < 2
        ME = MException('Matlab:WrongSize','Function requires more than one data point!');
        throw(ME);
    end
    if  ~exist('idx','var')
        if size(ps,2)>2 
            ME = MException('Matlab:VariableNotFound','No indices provided!');
            throw(ME);
        else %ps has exactly two objectives
            idx = [1,2];
        end
    else %idx variable provided
        if numel(idx) ~= 2 %how many indices are provided
                    ME = MException('Matlab:WrongSize','Function needs exactly two indices!');
            throw(ME);
        end

        if any(idx<1) || any(idx>size(ps,2)) %do the indices have the correct dimension?
                    ME = MException('Matlab:WrongIdx','Indices out of bounds!');
            throw(ME);
        end
    end

    
    figure
    scatter(ps(:,idx(1)),ps(:,idx(2)),'black')
    xlabel('f_1[a.u]')
    ylabel('f_2[a.u]')
end