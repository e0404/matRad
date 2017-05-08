function [ resultGUI ] = matRad_reCalcMCNRBExD(cst, dij,resultGUI )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


   dp       = resultGUI.physicalDose;
   ix       = dp > 0;
   LETd     = zeros(dij.numOfVoxels,1);
   RBEmax   = zeros(dij.numOfVoxels,1);
   RBEmin   = zeros(dij.numOfVoxels,1);
   LETd(ix) = (dij.mLETDose{1}(ix,:)  * resultGUI.w)./dp(ix);

   if ~isfield(dij,'alphaX')
      % Only take voxels inside patient.
      V = [cst{:,4}];
      V = unique(vertcat(V{:}));
      dij.alphaX = zeros(dij.numOfVoxels,1);
      dij.betaX  = zeros(dij.numOfVoxels,1);
      dij.abX    = zeros(dij.numOfVoxels,1);
      %set overlap priorities
      cst  = matRad_setOverlapPriorities(cst);

      for i = 1:size(cst,1)
        % find indices of structures related to V
        [~, row] = ismember(vertcat(cst{i,4}{:}),V,'rows'); 

        % check if cst is compatiable 
        if ~isempty(cst{i,5}) && isfield(cst{i,5},'alphaX') && isfield(cst{i,5},'betaX') 
            dij.alphaX(V(row)) = cst{i,5}.alphaX;
            dij.betaX(V(row))  = cst{i,5}.betaX;               
        else
            dij.alphaX(V(row)) = 0.1;    % default parameter  
            dij.betaX(V(row))  = 0.05;   % default parameter
            fprintf(['matRad: using default alpha_x and beta_x parameters for ' cst{i,2} ' \n']);
        end

      end
      dij.abX(dij.betaX>0) = dij.alphaX(dij.betaX>0)./dij.betaX(dij.betaX>0);
   end
   
   
   RBEmax(ix) = 0.999064 + ((0.35605  * LETd(ix) )./ dij.abX(ix));  
   RBEmin(ix) = 1.1012 + (-0.0038703  * real(sqrt(dij.abX(ix))) .* LETd(ix)); 
   resultGUI.RBExDose  = reshape(0.5 .* (sqrt(dij.abX.^2 + (4*dp(:).*dij.abX.*RBEmax) + (4*dp(:).^2 .* RBEmin.^2)) - dij.abX),dij.dimensions);


end

