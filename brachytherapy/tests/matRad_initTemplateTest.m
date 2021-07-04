classdef matRad_calcBrachyDoseTest < matlab.unittest.TestCase
    % tests all the functions called during dose calculation
     properties
        OriginalPath
    end
    
    methods (TestMethodSetup)
        function addBrachyToPath(testCase)
            testCase.OriginalPath = path;
            addpath(fileparts(fileparts(fileparts(mfilename('fullpath'))))); %matRad
            addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'phantoms')); %matRad/phantoms
        end
    end
    
    methods (TestMethodTeardown)
        function restorePath(testCase)
            path(testCase.OriginalPath);
        end
    end
    
    methods(Test)
        % test if rotation matrix yields the right result
        function rotMtx(testCase)
            HorPoints  = 8;
            VertPoints = 6;
            [templXmesh,templYmesh] = meshgrid(0:HorPoints-1,0:VertPoints-1);
            templX = reshape(templXmesh,[],1);
            templY = reshape(templYmesh,[],1);
            templZ = zeros(length(templX),1);
            stf.template.template2D = [templX';templY';templZ'];
            
            Xscale          = 20; % [mm]
            Yscale          = 20; % [mm]
            seedDistance      = 10; % [mm]
            seedsNo           = 11; 
            isoCenter                = matRad_getIsoCenter(cst,ct,0);
            %unit vectors of displaced, rotated template coordinate system
            xDir = normalize([1,0,0],'norm');
        	yDir = normalize([0,1,0],'norm');
            zDir = cross(pln.propStf.orientation.Xdir,pln.propStf.orientation.Ydir);
            offset = [-68,-47,-64]; 
            template3D = Xsc.*xDir'*temp2D(1,:) + Ysc.*yDir'*temp2D(2,:) + zDir'*temp2D(3,:)+ offs';
            stf.template.template3D = template3D;
            
            stf = matRad_generateBrachyStf(ct, cst, pln);
            dij = matRad_calcBrachyDose(ct,stf,pln,cst);
            testCase.verifyTrue(isfield(dij,'doseGrid'));
            testCase.verifyTrue(isfield(dij,'physicalDose'));
            testCase.verifyTrue(isfield(dij,'totalNumOfBixels'));
            testCase.verifyTrue(iscell(dij.physicalDose));
            

end




        template3D = Xsc.*xDir'*temp2D(1,:) + Ysc.*yDir'*temp2D(2,:) + zDir'*temp2D(3,:)+ offs';
        stf.template.template3D = template3D;