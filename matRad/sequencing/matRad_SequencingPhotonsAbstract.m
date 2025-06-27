classdef  (Abstract) matRad_SequencingPhotonsAbstract < matRad_SequencingBase

    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes her
    properties
        numOfLevels;
    end

    methods


        function [D_0,D_k, shapes,calFac, indInMx] = initBeam(this,stf, wCurr)

                numOfRaysPerBeam = stf.numOfRays; 
                X = ones(numOfRaysPerBeam,1)*NaN;
                Z = ones(numOfRaysPerBeam,1)*NaN;
                
                for j = 1:stf.numOfRays
                    X(j) = stf.ray(j).rayPos_bev(:,1);
                    Z(j) = stf.ray(j).rayPos_bev(:,3);
                end
                
                % sort bixels into matrix
                minX = min(X);
                maxX = max(X);
                minZ = min(Z);
                maxZ = max(Z);
                
                dimOfFluenceMxX = (maxX-minX)/stf.bixelWidth + 1;
                dimOfFluenceMxZ = (maxZ-minZ)/stf.bixelWidth + 1;

                %Create the fluence matrix.
                fluenceMx = zeros(dimOfFluenceMxZ,dimOfFluenceMxX);
                
                % Calculate X and Z position of every fluence's matrix spot z axis =
                % axis of leaf movement!
                xPos = (X-minX)/stf.bixelWidth+1;
                zPos = (Z-minZ)/stf.bixelWidth+1;
                
                % Make subscripts for fluence matrix
                indInMx = zPos + (xPos-1)*dimOfFluenceMxZ;
                
                %Save weights in fluence matrix.
                fluenceMx(indInMx) = wCurr;

                % Stratification
                calFac = max(fluenceMx(:));
                D_k = round(fluenceMx/calFac*this.numOfLevels); 
                
                % Save the stratification in the initial intensity matrix D_0.
                D_0 = D_k;

                % container to remember generated shapes; allocate space for 10000 shapes
                shapes = NaN*ones(dimOfFluenceMxZ,dimOfFluenceMxX,10000);
                
        end


    end
end