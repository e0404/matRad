%% Sparse_read_header
% Import a text file containing the dose at sparse format
%
%% Syntax
% |Sparse_info = Sparse_read_header(filename)|
%
%
%% Description
% |Sparse_info = Sparse_read_header(filename)| Description
%
%
%% Input arguments
% |filename| - _STRING_ - Name of the text file containing the dose information
%
%
%% Output arguments
%
% |Sparse_info| - _STRUCTURE_ - Data structure containing the sparse doose information
%
% * -- |Sparse_info.Header_file| - _STRING_ - Name of the file containing the dose information
% * -- |Sparse_info.Format| - _STRING_ - 'MCsquare Sparse Format'
% * -- |Sparse_info.NbrSpots| - _INTEGER_ - Number of spots
% * -- |Sparse_info.ImageSize| - _INTEGER_ - Number of pixel (x,y) in the image
% * -- |Sparse_info.VoxelSpacing| - _SCALAR_ - Dimension of the pixels (mm)
% * -- |Sparse_info.SimulationMode{1}| - __ -
%
%
%% Contributors
% Authors : K. Souris (open.reggui@gmail.com)

function Sparse_info = Sparse_read_header(filename)

  fid=fopen(filename,'r');
  if(fid < 0)
      error(['Unable to open file ' filename])
  end
  
  Sparse_info.Header_file = filename;
  Sparse_info.Format = 'MCsquare Sparse Format';
  
  while(~feof(fid))
      
      str=fgetl(fid);
      
      % find the index of '='
      s=find(str=='=',1,'first');
      
      if(~isempty(s))
        type=str(1:s-1);    % Copy the key
        data=str(s+1:end);  % Copy the value
        
        % Remove spaces
        while(type(end)==' '); type=type(1:end-1); end
        while(data(1)==' '); data=data(2:end); end
        
      else
        continue
      end
      
      switch(type)
          case 'NbrSpots'
              Sparse_info.NbrSpots = sscanf(data, '%d')';
          case 'ImageSize'
              Sparse_info.ImageSize = sscanf(data, '%d')';
          case 'VoxelSpacing'
              Sparse_info.VoxelSpacing = sscanf(data, '%f')';
          case 'SimulationMode'
              if(isfield(Sparse_info, 'SimulationMode'))
                Sparse_info.SimulationMode{numel(Sparse_info.SimulationMode)+1} = data;
              else
                Sparse_info.SimulationMode{1} = data;
              end
          otherwise
              Sparse_info.(type)=data;
      end

  end
  
  fclose(fid);
  
end
