function stf = matRad_bixelPhspVmc(stf,masterRayPosBEV,vmcOptions)

% after everything is done and working, add in code to verify if files
% already exist



switch vmcOptions.version
    case 'Carleton'
        phspPath = fullfile(fileparts(mfilename('fullpath')), 'run', 'phsp');
    case 'dkfz'
        phspPath = fullfile(fileparts(mfilename('fullpath')), 'runs', 'phsp');
end

fname_full  = fullfile(phspPath,sprintf('%s.egsphsp1',vmcOptions.phspBaseName));
fid_full    = fopen(fname_full,'r');

SAD2SCD = vmcOptions.SCD./stf(1).SAD;
% in cm
bixelWidth  = stf(1).bixelWidth.*SAD2SCD/10;
X           = masterRayPosBEV(:,1).*SAD2SCD/10;
Y           = -masterRayPosBEV(:,3).*SAD2SCD/10; % minus sign necessary since to get from BEAM coord. to DICOM, we do a rotation, NOT reflection

% determine file names, check for existence
numBixels       = size(masterRayPosBEV,1);
fname_bixels    = cell(numBixels,1);
writeFiles      = false;
for i = 1:numBixels
    
    % file name
    fname_bixels{i}         = fullfile(phspPath,sprintf('%s_bixelWidth%f_X%f_Y%f.egsphsp1',vmcOptions.phspBaseName,bixelWidth,X(i),Y(i)));
    
    if ~exist(fname_bixels{i},'file')
        % if any file doesn't exist, then we want to write new phsp file
        writeFiles = true;
    end
end

% give file name to ray
for i = 1:numel(stf)
    
    for j = 1:stf(i).numOfRays
        
        % find correct bixel
        bixelInd = all(repelem(stf(i).ray(j).rayPos_bev,numBixels,1) == masterRayPosBEV,2);
        % write filename
        stf(i).ray(j).phspFileName = fname_bixels{bixelInd};
    end
end

% FIX THIS TO GENERATE PHSP FILES FOR ALL BIXELS IN FIELD

if writeFiles
    % only do read/write files if they don't already exist
    
    %% extract information from full phsp file, write to bixel files
    
    % set up arrays for the bixel phsp files
    fid_bixels          = cell(numBixels,1);
    header_bixels       = cell(numBixels,1);
    firstParticle_bixels   = false(numBixels,1);
    
    % open header of full phsp
    [fid_full, header_full] = getHeader(fid_full);
    mode = char(header_full.MODE_RW(5));
    
    % loop through each record in full phsp
    fprintf('matRad: creating bixel phsp files... ');
    for i = 1:header_full.NPPHSP
        
        % extract record
        [fid_full, record] = getRecord(fid_full,mode);
        
        % sort into correct bixel
        bixelInd = find(sum(abs([X Y]-repelem([record.X record.Y],numBixels,1)) < repelem(bixelWidth/2,numBixels,2),2) == 2,1,'first');
        
        if ~isempty(bixelInd)
            
            if isempty(fid_bixels{bixelInd})
                % if not previously opened, open the file, write the header
                % the header is just a dummy for now
                header_bixels{bixelInd}             = header_full;
                % these variables will change throughout the read/write process
                header_bixels{bixelInd}.NPPHSP      = 0;
                header_bixels{bixelInd}.NPHOTPHSP   = 0;
                header_bixels{bixelInd}.EKMAXPHSP   = 0;
                header_bixels{bixelInd}.EKMINPHSPE  = 1000;
                
                % open file, write header
                fid_bixels{bixelInd} = fopen(fname_bixels{bixelInd},'W'); % turn this to 'W'?
                writeHeader(fid_bixels{bixelInd},header_bixels{bixelInd});
            end
            
            %% bixel header
            
            % increment number of particles
            header_bixels{bixelInd}.NPPHSP = header_bixels{bixelInd}.NPPHSP+1;
            
            % modify max/min energies, increment number of photons
            % must determine particle type using LATCH
            LATCH = de2bi(record.LATCH,32);
            if LATCH(30:31) == [0 0]
                % photon
                header_bixels{bixelInd}.EKMAXPHSP   = max(header_bixels{bixelInd}.EKMAXPHSP,abs(record.E));
                header_bixels{bixelInd}.NPHOTPHSP   = header_bixels{bixelInd}.NPHOTPHSP+1;
                
            elseif LATCH(30:31) == [0 1]
                % electron
                header_bixels{bixelInd}.EKMINPHSPE  = min(header_bixels{bixelInd}.EKMINPHSPE,abs(record.E)-0.511);
                header_bixels{bixelInd}.EKMAXPHSP   = max(header_bixels{bixelInd}.EKMAXPHSP,abs(record.E)-0.511);
                
            elseif LATCH(30:31) == [1 0]
                % positron
                header_bixels{bixelInd}.EKMINPHSPE  = min(header_bixels{bixelInd}.EKMINPHSPE,abs(record.E)-0.511);
                header_bixels{bixelInd}.EKMAXPHSP   = max(header_bixels{bixelInd}.EKMAXPHSP,abs(record.E)-0.511);
                
            elseif LATCH(30:31) == [1 1]
                error('Electron and positron???')
                
            end
            
            %% bixel record
            
            % is this the first particle scored from a new primary history?
            if record.E < 0
                % if it is, then we want ALL bixel phsp files to reflect this
                % i.e., the next particle in all bixel phsp files should have
                % a negative energy (first particle scored from a new primary
                % history)
                
                firstParticle_bixels(:) = true;
            end
            
            if firstParticle_bixels(bixelInd)
                % this is the first from a new primary history
                % make the energy negative
                record.E = -abs(record.E);
                % make sure next particle is not negative
                firstParticle_bixels(bixelInd) = false;
            else
                % these particles should already have positive energy
                if record.E < 0
                    % SHOULDN'T HAPPEN
                    warning('NEGATIVE ENERGY')
                end
            end
            
            % write the record
            writeRecord(fid_bixels{bixelInd},record,mode);
        end
        
        % display progress
        if mod(i,max(1,round(header_full.NPPHSP/200))) == 0
            matRad_progress(i/max(1,round(header_full.NPPHSP/200)),...
                floor(header_full.NPPHSP/max(1,round(header_full.NPPHSP/200))));
        end
    end
    
    
    %% clean up bixel phsp headers
    fprintf('matRad: updating bixel phsp file headers... ');
    for i = 1:numBixels
        
        if header_bixels{i}.EKMINPHSPE == 1000
            header_bixels{i}.EKMINPHSPE = 0;
        end
        
        % seek to beginning of file
        fseek(fid_bixels{i},0,'bof');
        
        % write updated header
        writeHeader(fid_bixels{i},header_bixels{i});
        
        % close file
        fclose(fid_bixels{i});
    end
    fprintf('Done!\n')
    
end

end


% read/write functions

function [fid, header] = getHeader(fid)

header.MODE_RW      = fread(fid,5,'uint8');
header.NPPHSP       = fread(fid,1,'int32');
header.NPHOTPHSP    = fread(fid,1,'int32');
header.EKMAXPHSP    = fread(fid,1,'float32');
header.EKMINPHSPE   = fread(fid,1,'float32');
header.NINCPHSP     = fread(fid,1,'float32');
header.garbage      = fread(fid,3,'int8');

end

function [fid, record] = getRecord(fid,mode)

record.LATCH    = fread(fid,1,'uint32');
record.E        = fread(fid,1,'float32');
record.X        = fread(fid,1,'float32');
record.Y        = fread(fid,1,'float32');
record.U        = fread(fid,1,'float32');
record.V        = fread(fid,1,'float32');
record.WT       = fread(fid,1,'float32');

if mode == 2
    record.ZLAST    = fread(fid,1,'float32');
end

end

function writeHeader(fid,header)

fwrite(fid,header.MODE_RW,'uint8');
fwrite(fid,header.NPPHSP,'int32');
fwrite(fid,header.NPHOTPHSP,'int32');
fwrite(fid,header.EKMAXPHSP,'float32');
fwrite(fid,header.EKMINPHSPE,'float32');
fwrite(fid,header.NINCPHSP,'float32');
fwrite(fid,header.garbage,'int8');

end

function writeRecord(fid,record,mode)

fwrite(fid,record.LATCH,'uint32');
fwrite(fid,record.E,'float32');
fwrite(fid,record.X,'float32');
fwrite(fid,record.Y,'float32');
fwrite(fid,record.U,'float32');
fwrite(fid,record.V,'float32');
fwrite(fid,record.WT,'float32');

if mode == 2
    fwrite(fid,record.ZLAST,'float32');
end

end