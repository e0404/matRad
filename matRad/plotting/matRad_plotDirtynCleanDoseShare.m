function [stackedbarDose,physDoseInDepth,RBExDcurve,RBEcurve,h,d] = matRad_plotDirtynCleanDoseShare(definedEnd,index,ct,ctCube,dij,resultGUI,add,LET_thres,k,displayComparison)
%
% call
%   [LETbeamletSpectrum, PhysDose_LET] = matRad_DoseComparison(definedEnd,index,ct,ctCube,dij,resultGUI,add,LET_thres,k)
%   [LETbeamletSpectrum, PhysDose_LET] = matRad_DoseComparison(definedEnd,index,ct,ctCube,dij,resultGUI,add,LET_thres,k,displayComparison)
%
% input
%   definedEnd:             last voxel to consider
%   index:                  voxel coordinates from the cube index in matRadGUI in [y x z]
%   ct:                     matRad ct struct
%   ctCube:                 matRad ct.cube in ct struct
%   dij:                    matRad dij struct
%   resultGUI:              matRad resultGUI struct
%   add:                    add [y x z] to your index: 0 for no change, 1 for one step and so on
%   LET_thres:              LET threshold: above this = dirty dose, below this = clean dose
%   k:                      k the dimension that should change
%   displayComparisony:     displays a comparison of dose and RBE
%
% output
%   stackedbarDose:         bar plot with stacked dose shares adding up to the total dose
%   physDoseInDepth:        total physical dose plotted in a line depending on the penetration depth
%   RBExDcurve:             curve to show the RBExD distribution of different voxels
%   RBEcurve:               curve to show the RBE distribution of different voxels
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initialize arrays
    h = [];
    l = [];
    t = [];
    d = [];
    RBE = [];
    RBExD = [];
    
    for i = 1:definedEnd                                                % loop to describe a line of voxels within a certain area defined by the index
        
        % calculating dirty and clean dose for each voxel; don't display plots of matRad_returnDirtyandCleanDose, use []
        [highLETphysDose,lowLETphysDose,totalphysDose] = matRad_returnDirtyandCleanDose(index,ct,ctCube,dij,resultGUI,LET_thres,[],[],[]);
        
        depth = index(k);
      
        % safe results
        h = [h highLETphysDose];
        l = [l lowLETphysDose];
        t = [t totalphysDose];
        d = [d depth];
        
        % calculating RBExD and RBE for each voxel and storing them into arrays
        rowRBExD = resultGUI.RBExD(index(1),index(2),index(3));
        rowRBE = resultGUI.RBE(index(1),index(2),index(3));

        RBExD = [RBExD rowRBExD];
        RBE = [RBE rowRBE];
        
        index = index + add;                                            % describes next voxel
    end     

    if ~exist('displayComparison','var') || isempty(displayComparison)  % for showing image, displayComparison must be 1
        displayComparison = [];
    elseif displayComparison == 1
            sumDose = [h' l'];                                          % stores the dirty and clean dose into one array
            figure;
            hold on
            stackedbarDose = bar(d,sumDose,"stacked",'EdgeColor',[0 1 0],'LineWidth',1);    % plots stacked bars of dirty and clean dose to see the respective share of the total dose
            physDoseInDepth = plot(d,t,"Color",[0 0 1],'LineWidth',1);                      % plots the total physical dose
            RBExDcurve = plot(d,RBExD, "Color", [1 0 0], 'LineWidth', 1);                   % plots the RBExD
            RBEcurve = plot(d,RBE, "Color", [1 0 1], 'LineWidth', 1);                       % plots the RBE
            xlabel('depth in mm'); ylabel('physical dose in Gy');
            legend('dirty share of physical dose','clean share of physical dose','total physical dose', 'RBExD','RBE',Location='best');
    end

end