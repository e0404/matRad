function cube = matRad_exportCtTOPAS(ct, runsPath, basematerial)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad export CT RSP data for TOPAS simulation
%
% call
%   matRad_exportCtTOPAS(ct, path, material)
%
% input
%   ct:             ct cube
%   runsPath:       path where to save the files for MC simulation
%   basematerial:   base material to be scaled to corresponding RSP
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%arrayOrdering='C'; % row-major order
arrayOrdering='F'; % column-major order
% arrayOrdering='Fotran' % column-major order
vecRSPlength=10000;

% the image cube contains the indexing of materials
% since at the moment TOPAS does not support ushort
% the materials should have indexes between 0 and 32767
% therefore, the maximum length of the vector is 32768
%answer = ''
%while isempty(answer)
%  prompt = 'Length of the imaging vector [2-32768]:';
%  answer = input(prompt)
%  if answer > 1 && answer <= 32768
%    vecRSPlength=answer
%  else
%    answer = ''
%  end
%end

medium = basematerial;

if isequal(arrayOrdering,'C')
    disp('Exporting cube in C ordering...')
    permutation = [3 1 2];
else
    disp('Exporting cube in FORTRAN ordering...')
    permutation = [2 1 3];
end

cubeExport = 'HUSchneiderToWater'; %'RSP'; %'HUSchneiderToWater';
cubeExport = 'RSP';
checkMaterial = false;
rspCubeMethod = 2; %1: TsBox with variable voxels, 2: TsImageCube with Density Bins and Custom Image converter


paramFile = 'matRad_cube.txt';
dataFile = 'matRad_cube.dat';

nbVoxels = prod(ct.cubeDim);

outfile = fullfile(runsPath, paramFile);
disp(['Writing data to ',outfile])
h = fopen(outfile,'w+');

switch cubeExport
    case 'RSP'
        rspCube = ct.cube{1};
        rspCube = permute(rspCube,permutation); %  X,Y,Z ordering
                      
        fbase = fopen(['materials/' medium '.txt'],'r');
        while ~feof(fbase)
            strLine = fgets(fbase); %# read line by line
            fprintf(h,'%s',strLine);
        end
        fclose(fbase);
        
        minRSP = min(rspCube(:));
        maxRSP = max(rspCube(:));
        
        % to avoid zero density
        if minRSP<1.e-6
            minRSP=1.e-6;
        end
        
        % case of homogenous water medium (i.e., RSP=1)
        if (length(unique(ct.cube{1})) == 1) && (unique(ct.cube{1} == 1))
            fprintf(h,'s:Ge/Patient/Parent="World"\n');
            fprintf(h,'s:Ge/Patient/Type= "TsBox"\n');
            fprintf(h,'s:Ge/Patient/Material = "%s"\n',medium);
            fprintf(h,'d:Ge/Patient/HLX      = %f mm\n',0.5*ct.cubeDim(2)*ct.resolution.x);
            fprintf(h,'d:Ge/Patient/HLY      = %f mm\n',0.5*ct.cubeDim(1)*ct.resolution.y);
            fprintf(h,'d:Ge/Patient/HLZ      = %f mm\n',0.5*ct.cubeDim(3)*ct.resolution.z);
            fprintf(h,'i:Ge/Patient/XBins    = %d\n',ct.cubeDim(2));
            fprintf(h,'i:Ge/Patient/YBins    = %d\n',ct.cubeDim(1));
            fprintf(h,'i:Ge/Patient/ZBins    = %d\n',ct.cubeDim(3));
            
            % otherwise
        else
            
            % to avoid issues with homogenous media
            if minRSP+1.e-6>maxRSP
                warning('Use only one RSP value')
                vecRSPlength = 2;
                minRSP = 0.5*maxRSP;
            end
            
            dRSP = (maxRSP-minRSP)/(vecRSPlength-1);
            upperRSP = maxRSP+dRSP;
            
            ixRSP = round((rspCube-minRSP)/dRSP)+1;
            
            fprintf(h,'s:Ge/Patient/Parent="World"\n');
            fprintf(h,'i:Ma/%s/VariableDensityBins = %d\n',medium,vecRSPlength);
            fprintf(h,'u:Ma/%s/VariableDensityMin = %f\n',medium,minRSP);
            fprintf(h,'u:Ma/%s/VariableDensityMax = %f\n',medium,upperRSP);
            
            if rspCubeMethod == 1
                fprintf(h,'s:Ge/Patient/Type= "TsBox"\n');
                fprintf(h,'s:Ge/Patient/Material = "%s"\n',medium);
                fprintf(h,'d:Ge/Patient/HLX      = %f mm\n',0.5*ct.cubeDim(2)*ct.resolution.x);
                fprintf(h,'d:Ge/Patient/HLY      = %f mm\n',0.5*ct.cubeDim(1)*ct.resolution.y);
                fprintf(h,'d:Ge/Patient/HLZ      = %f mm\n',0.5*ct.cubeDim(3)*ct.resolution.z);
                fprintf(h,'i:Ge/Patient/XBins    = %d\n',ct.cubeDim(2));
                fprintf(h,'i:Ge/Patient/YBins    = %d\n',ct.cubeDim(1));
                fprintf(h,'i:Ge/Patient/ZBins    = %d\n',ct.cubeDim(3));
                fprintf(h,'sv:Ge/Patient/VoxelMaterials = %d\n',nbVoxels);
                
                voxelString = num2str(ixRSP(:)'-1,['"' medium '_VariableDensityBin_%d"\n']);
                
                %for ix=1:nbVoxels
                %    fprintf(h,'"%s"\n',[ medium '_VariableDensityBin_' num2str(ixRSP(ix)-1)]);
                %end
                fprintf(h,voxelString);
                
                if checkMaterial
                    for ix=1:nbVoxels
                        rspMaterial{ix} = [ medium '_VariableDensityBin_' num2str(ixRSP(ix)-1)];
                    end
                    materialsUsed = unique(rspMaterial);
                    
                    %      fprintf(h,'sv:Sc/ExtractData/Material = %d ',length(materialsUsed))
                    %      for ix=1:length(materialsUsed)
                    %      fprintf(h,'"%s" ',materialsUsed{ix});
                    %      end
                    %      fprintf(h,'\n')
                end
                fclose(h);
                cube = rspCube;
                
            elseif rspCubeMethod == 2
                fprintf(h,'s:Ge/Patient/Type = "TsImageCube"\n');
                fprintf(h,'b:Ge/Patient/DumpImagingValues = "True"\n');
                fprintf(h,'s:Ge/Patient/BaseMaterial = "%s"\n',medium);
                fprintf(h,'i:Ge/Patient/MaterialIxMax = %d\n',vecRSPlength);
                fprintf(h,'s:Ge/Patient/InputDirectory = "./"\n');
                fprintf(h,'s:Ge/Patient/InputFile = "%s"\n',dataFile);
                fprintf(h,'s:Ge/Patient/ImagingtoMaterialConverter = "matrad"\n');
                fprintf(h,'i:Ge/Patient/NumberOfVoxelsX = %d\n',ct.cubeDim(2));
                fprintf(h,'i:Ge/Patient/NumberOfVoxelsY = %d\n',ct.cubeDim(1));
                fprintf(h,'iv:Ge/Patient/NumberOfVoxelsZ = 1 %d\n',ct.cubeDim(3));
                fprintf(h,'d:Ge/Patient/VoxelSizeX       = %.3f mm\n',ct.resolution.x);
                fprintf(h,'d:Ge/Patient/VoxelSizeY       = %.3f mm\n',ct.resolution.y);
                fprintf(h,'dv:Ge/Patient/VoxelSizeZ       = 1 %.3f mm\n',ct.resolution.z);
                fprintf(h,'s:Ge/Patient/DataType  = "SHORT"\n');
                fclose(h);
                
                % write data
                h = fopen(fullfile(runsPath, dataFile),'w');
                fwrite(h,ixRSP(:)-1,'short');
                fclose(h);
                cube = rspCube;
                
                %Water equivalent Schneider Converter
            elseif methodCube == 3

            end
        end
        
    case 'HUSchneiderToWater'        
        huCube = int32(permute(ct.cubeHU{1},permutation));
        
        rspHlut = matRad_readHLUT('matRad_default.hlut');
        densityCorrection = [];
        for i = 1:size(rspHlut,1)-1
            startVal = rspHlut(i,1);
            endVal = rspHlut(i+1,1);
            range = startVal:1:endVal-1;
            densityCorrection(end+1:end+numel(range)) = matRad_interp1(rspHlut(:,1),rspHlut(:,2),range);
        end
        densityCorrection(end+1) = rspHlut(end,2); %add last missing value
        
        
        %Write the Schneider Converter
        fprintf(h,'dv:Ge/Patient/DensityCorrection = %d %s %s\n',numel(densityCorrection),num2str(densityCorrection,'%f '),'g/cm3');
        fprintf(h,'iv:Ge/Patient/SchneiderHounsfieldUnitSections = 2 %d %d\n',rspHlut(1,1),rspHlut(end,1)+1);
        fprintf(h,'uv:Ge/Patient/SchneiderDensityOffset = 1 1\n');
        fprintf(h,'uv:Ge/Patient/SchneiderDensityFactor = 1 0\n');
        fprintf(h,'uv:Ge/Patient/SchneiderDensityFactorOffset = 1 %d\n\n',-rspHlut(1,1));
        fprintf(h,'i:Ge/Patient/MinImagingValue = %d\n',rspHlut(1,1));
        
        fprintf(h,'iv:Ge/Patient/SchneiderHUToMaterialSections = 2 %d %d\n',rspHlut(1,1),rspHlut(end,1)+1);
        fprintf(h,'sv:Ge/Patient/SchneiderElements = 2 "Hydrogen" "Oxygen"\n');
        fprintf(h,'uv:Ge/Patient/SchneiderMaterialsWeight1 = 2 0.111894 0.888106\n');
        
        %Write the Patient
        fprintf(h,'s:Ge/Patient/Parent="World"\n');
        fprintf(h,'s:Ge/Patient/Type = "TsImageCube"\n');
        fprintf(h,'b:Ge/Patient/DumpImagingValues = "True"\n');
        %fprintf(h,'s:Ge/Patient/BaseMaterial = "%s"\n',medium);
        %fprintf(h,'i:Ge/Patient/MaterialIxMax = %d\n',vecRSPlength);
        fprintf(h,'s:Ge/Patient/InputDirectory = "./"\n');
        fprintf(h,'s:Ge/Patient/InputFile = "%s"\n',dataFile);
        fprintf(h,'s:Ge/Patient/ImagingtoMaterialConverter = "Schneider"\n');
        fprintf(h,'i:Ge/Patient/NumberOfVoxelsX = %d\n',ct.cubeDim(2));
        fprintf(h,'i:Ge/Patient/NumberOfVoxelsY = %d\n',ct.cubeDim(1));
        fprintf(h,'iv:Ge/Patient/NumberOfVoxelsZ = 1 %d\n',ct.cubeDim(3));
        fprintf(h,'d:Ge/Patient/VoxelSizeX       = %.3f mm\n',ct.resolution.x);
        fprintf(h,'d:Ge/Patient/VoxelSizeY       = %.3f mm\n',ct.resolution.y);
        fprintf(h,'dv:Ge/Patient/VoxelSizeZ       = 1 %.3f mm\n',ct.resolution.z);
        fprintf(h,'s:Ge/Patient/DataType  = "SHORT"\n');
        %fprintf(h,'includeFile = HUtoMaterialSchneiderWater.txt');
        fclose(h);
        
        % write data
        h = fopen(fullfile(runsPath, dataFile),'w');
        fwrite(h,huCube,'short');
        fclose(h);
        cube = huCube;
    otherwise 
        disp('Cube Export not defined!');
end
end
