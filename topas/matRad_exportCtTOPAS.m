function [rspCube] = matRad_exportCtTOPAS(ct, runsPath, basematerial)
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
methodCube = 2;
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
exportRSPcube = true;
checkMaterial = false;

if exportRSPcube

  rspCube = ct.cube{1};
  rspCube = permute(rspCube,[2 1 3]); %  X,Y,Z ordering
  nbVoxels = ct.cubeDim(2)*ct.cubeDim(1)*ct.cubeDim(3);

%    rspCube = uint16(rspCube*1000.);

  if isequal(arrayOrdering,'C')
    disp('Exporting RSP cube in C ordering...')
    rspCube = permute(rspCube,[3 2 1]);
  else
    disp('Exporting RSP cube in FORTRAN ordering...')
  end

  outfile = fullfile(runsPath, 'matRad_RSPcube.txt');
  disp(['Writing data to ',outfile])
  h = fopen(outfile,'w+');

  warning(['Using base material as ' medium])
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

    if methodCube == 1
      fprintf(h,'s:Ge/Patient/Type= "TsBox"\n');
      fprintf(h,'s:Ge/Patient/Material = "%s"\n',medium);
      fprintf(h,'d:Ge/Patient/HLX      = %f mm\n',0.5*ct.cubeDim(2)*ct.resolution.x);
      fprintf(h,'d:Ge/Patient/HLY      = %f mm\n',0.5*ct.cubeDim(1)*ct.resolution.y);
      fprintf(h,'d:Ge/Patient/HLZ      = %f mm\n',0.5*ct.cubeDim(3)*ct.resolution.z);
      fprintf(h,'i:Ge/Patient/XBins    = %d\n',ct.cubeDim(2));
      fprintf(h,'i:Ge/Patient/YBins    = %d\n',ct.cubeDim(1));
      fprintf(h,'i:Ge/Patient/ZBins    = %d\n',ct.cubeDim(3));
      fprintf(h,'sv:Ge/Patient/VoxelMaterials = %d\n',nbVoxels);
      for ix=1:nbVoxels
        fprintf(h,'"%s"\n',[ medium '_VariableDensityBin_' num2str(ixRSP(ix)-1)]);
      end

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

    elseif methodCube == 2
      fprintf(h,'s:Ge/Patient/Type = "TsImageCube"\n');
      fprintf(h,'b:Ge/Patient/DumpImagingValues = "True"\n');
      fprintf(h,'s:Ge/Patient/BaseMaterial = "%s"\n',medium);
      fprintf(h,'i:Ge/Patient/MaterialIxMax = %d\n',vecRSPlength);
      fprintf(h,'s:Ge/Patient/InputDirectory = "./"\n');
      fprintf(h,'s:Ge/Patient/InputFile = "%s"\n','matRad_RSPcube.dat'));
      fprintf(h,'s:Ge/Patient/ImagingtoMaterialConverter = "matrad"\n');
      fprintf(h,'i:Ge/Patient/NumberOfVoxelsX = %d\n',ct.cubeDim(2));
      fprintf(h,'i:Ge/Patient/NumberOfVoxelsY = %d\n',ct.cubeDim(1));
      fprintf(h,'i:Ge/Patient/NumberOfVoxelsZ = %d\n',ct.cubeDim(3));
      fprintf(h,'d:Ge/Patient/VoxelSizeX       = %.3f mm\n',ct.resolution.x);
      fprintf(h,'d:Ge/Patient/VoxelSizeY       = %.3f mm\n',ct.resolution.y);
      fprintf(h,'d:Ge/Patient/VoxelSizeZ       = %.3f mm\n',ct.resolution.z);
      fprintf(h,'s:Ge/Patient/DataType  = "SHORT"\n');
      fclose(h);

      % write data
      h = fopen(fullfile(runsPath, 'matRad_RSPcube.dat'),'w');
      fwrite(h,ixRSP(:)-1,'short');
      fclose(h);
    end
  end
end
