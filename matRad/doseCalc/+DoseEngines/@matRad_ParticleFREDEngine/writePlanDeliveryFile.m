function writePlanDeliveryFile(this, fName)
% FRED helper to write file for plan delivery routine
% call
%   writePlanDeliveryFile(fName)
% 
% input
%   fName: tring specifying the file path and name for saving the data.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2023 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matRad_cfg = MatRad_Config.instance();

fID = fopen(fName, 'w');
try
    fprintf(fID, '#Include file defining fields and layers geometry\n');
    fprintf(fID, 'include: inp/plan/plan.inp\n');
    fprintf(fID, '\n');
    

    fprintf(fID, '#Define the fields\n');
    fprintf(fID, 'for(currField in plan.get(''Fields''))<\n');
        fprintf(fID, '\tfield<\n');
            fprintf(fID, '\t\tID = ${currField.get(''fieldNumber'')}\n');
            fprintf(fID, '\t\tO = [0,${plan.get(''SAD'')},0]\n');
            fprintf(fID, '\t\tL = ${currField.get(''dim'')}\n');
            fprintf(fID, '\t\tpivot = [0.5,0.5,0.5]\n');

            fprintf(fID, '\t\tl = [0, 0, -1]\n');
            fprintf(fID, '\t\tu = [1, 0 ,0]\n');
           
        fprintf(fID, '\tfield>\n');

        fprintf(fID, '\n');
        fprintf(fID, '\t#Deactivate the fields to avoid geometrical overlap\n');
        fprintf(fID, '\tdeactivate: field_${currField.get(''fieldNumber'')}\n');
    fprintf(fID, 'for>\n\n');

    %loop over fields
    fprintf(fID, 'for(currField in plan.get(''Fields''))<\n');
        fprintf(fID,'\n');
        fprintf(fID, '\tdef: fieldIdx = currField.get(''fieldNumber'')\n');
        fprintf(fID,'\n');

        % activate current filed            
        fprintf(fID, '\t#Activate current field\n');
        fprintf(fID, '\tactivate: field_$fieldIdx\n');
        fprintf(fID,'\n');

        % Gantry/COunch angles
        fprintf(fID,'\t#Collect Gantry and Couch angles\n');
        fprintf(fID, '\tdef: GA = currField.get(''GA'')\n');
        fprintf(fID, '\tdef: CA = currField.get(''CA'')\n');
        fprintf(fID,'\n');

        % Isocenter
        fprintf(fID, '\t#Collect Isocenter\n');
        fprintf(fID, '\tdef: ISO = currField.get(''ISO'')\n');
        fprintf(fID, '\n');

        fprintf(fID, '\t#First move the patient so that the Isocenter is now in the center of the Room coordinate system\n');
        fprintf(fID, '\ttransform: Phantom move_to ${ISO.item(0)} ${ISO.item(1)} ${ISO.item(2)} Room\n');
        fprintf(fID, '\n');

        fprintf(fID, '\t#Second rotate the patient according to the gantry and couch angles.\n');
        fprintf(fID, '\t#In this configuration the fileds are always fixed in +SAD in y direction and the patient is rotated accordingly\n');
        fprintf(fID, '\ttransform: Phantom rotate y ${CA} Room\n');
        fprintf(fID, '\ttransform: Phantom rotate z ${GA} Room\n');
        fprintf(fID, '\n');

        fprintf(fID, '\tfor(layer in currField.get(''Layers''))<\n');
        fprintf(fID, '\n');

            fprintf(fID, '\t\t#Recover parameters of the current energy layer\n');
            fprintf(fID, '\t\tdef: currEnergy  = layer.get(''Energy'')\n');
            fprintf(fID, '\t\tdef: currEspread = layer.get(''Espread'')\n');

            switch this.sourceModel
                case 'gaussian'
                    fprintf(fID, '\t\tdef: currFWHM   = layer.get(''FWHM'')\n');

                case 'emittance'
                    fprintf(fID, '\t\tdef: currEmittanceX = layer.get(''emittanceX'')\n');

                    fprintf(fID, '\t\tdef: currTwissAlphaX = layer.get(''twissAlphaX'')\n');
                    fprintf(fID, '\t\tdef: currTwissBetaX = layer.get(''twissBetaX'')\n');
                    fprintf(fID, '\t\tdef: currReferencePlane = layer.get(''emittanceRefPlaneDistance'')\n');
                case 'sigmaSqrModel'
                    fprintf(fID, '\t\tdef: currSQr_a = layer.get(''sSQr_a'')\n');
                    fprintf(fID, '\t\tdef: currSQr_b = layer.get(''sSQr_b'')\n');
                    fprintf(fID, '\t\tdef: currSQr_c = layer.get(''sSQr_c'')\n');
            end

            fprintf(fID, '\n');
            fprintf(fID, '\t\tfor(beamlet in layer.get(''beamlets''))<\n');
                fprintf(fID, '\t\t\tpb<\n');
                    fprintf(fID, '\t\t\t\tID      = ${beamlet.get(''beamletID'')}\n');
                    fprintf(fID, '\t\t\t\tfieldID = $fieldIdx\n');
                    switch this.machine.meta.radiationMode
                        case 'protons'
                            fprintf(fID, '\t\t\t\tparticle = proton\n');
                        case 'carbon'
                            fprintf(fID, '\t\t\t\tparticle = C12\n');
                    end
                    fprintf(fID, '\t\t\t\tT	= $currEnergy\n');
                    fprintf(fID, '\t\t\t\tEFWHM	= $currEspread\n');
        
                    switch this.sourceModel
                        case 'gaussian'
        
                            fprintf(fID, '\t\t\t\tXsec = gauss\n');
                            fprintf(fID, '\t\t\t\tFWHM 	= $currFWHM\n');
                        case 'emittance'
                            fprintf(fID, '\t\t\t\tXsec = emittance\n');
                            fprintf(fID, '\t\t\t\temittanceX  = $currEmittanceX\n');
                            fprintf(fID, '\t\t\t\ttwissAlphaX = $currTwissAlphaX\n');
                            fprintf(fID, '\t\t\t\ttwissBetaX  = $currTwissBetaX\n');
                            fprintf(fID, '\t\t\t\temittanceRefPlaneDistance = 100\n');
        
                        case 'sigmaSqrModel'
                            fprintf(fID, '\t\t\t\tXsec = emittance\n');
                            fprintf(fID, '\t\t\t\tsigmaSqrModel = [${plan.get(''SAD'')},${currSQr_a},${currSQr_b}, ${currSQr_c}]\n');
                    end
        
                    fprintf(fID, '\n');
                    fprintf(fID, '\t\t\t\tP 	= ${beamlet.get(''P'')}\n');
                    fprintf(fID, '\t\t\t\tv 	= ${beamlet.get(''v'')}\n');
                    fprintf(fID, '\t\t\t\tN 	= ${beamlet.get(''w'')}\n');
                fprintf(fID, '\t\t\tpb>\n');
            fprintf(fID, '\t\tfor>\n');
        fprintf(fID, '\tfor>\n');
        fprintf(fID, '\n');
       
        fprintf(fID, '\t#Deliver all the pecil beams in this field\n');
        fprintf(fID, '\tdeliver: field_$fieldIdx\n');
        fprintf(fID, '\n');
         
        fprintf(fID, '\t#Deactivate the current field\n');
        fprintf(fID, '\tdeactivate: field_$fieldIdx\n');
        fprintf(fID, '\n');

        fprintf(fID, '\t#Restore the patient to original position\n');
        fprintf(fID, '\ttransform: Phantom rotate z ${-1*GA} Room\n');
        fprintf(fID, '\ttransform: Phantom rotate y ${-1*CA} Room\n');
        fprintf(fID, '\ttransform: Phantom move_to 0 0 0 Room\n');

        fprintf(fID, 'for>\n\n');
catch
    matRad_cfg.dispError('Failed to write planDelivery file');
end

fclose(fID);
end