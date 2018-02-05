function matRad_visSpotWeights(stf,weights)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualise spot weights per energy slice (or fluence map for photons respectively)
%   for single beams
% 
% call
%    matRad_visSpotWeights(stf,resultGUI.w)
%
% input
%   stf:              matRad stf struct
%   weights:          spot weights for bixels (resultGUI.w)

% output 
%      plots for all beams - for particles: scroll mousewheel to access different energies  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

numOfBeams = size(stf,2);

% construct look up table weight --> bixel
counter = 0;
for i = 1:numOfBeams
    for j = 1:stf(i).numOfRays
        for k = 1:stf(i).numOfBixelsPerRay(j)
            counter = counter + 1;
            bixelLut.beamNum(counter) = i;
            bixelLut.rayNum(counter) = j;
            bixelLut.bixelNum(counter) = k;
            bixelLut.energy(counter) = stf(i).ray(j).energy(k);
        end 
    end
end

    if  strcmp(stf(1).radiationMode, 'photons')
        
        % plot weights for all beams in one figure
        figure;
        rowNum = floor(sqrt(numOfBeams)); % number of rows in subplot
        colNum = ceil(numOfBeams/rowNum); % numberof columns in subplot
        % maximum weight
        w_max = max(weights);
        
        for i=1:numOfBeams
           
           % find range of ray positions within beam
           rayPos_mat = vertcat(stf(i).ray(:).rayPos_bev);
           x_min = min(rayPos_mat(:, 1));
           z_min = min(rayPos_mat(:, 3));
           x_max = max(rayPos_mat(:, 1));
           z_max = max(rayPos_mat(:, 3));

           % initialise weight matrix
           weight_matrix{i} = zeros((x_max-x_min)/stf(i).bixelWidth,(z_max-z_min)/stf(i).bixelWidth);

           for j=1:stf(i).numOfRays
               
               % find weight for this ray
               weight_ix = and((bixelLut.beamNum==i), (bixelLut.rayNum==j));
               % normalize weight
               w_norm = weights(weight_ix)/w_max;
               % find position in matrix corresponding to ray
               xpos = (stf(i).ray(j).rayPos_bev(1) + abs(x_min))/stf(i).bixelWidth +1;
               zpos = (stf(i).ray(j).rayPos_bev(3) + abs(z_min))/stf(i).bixelWidth +1;
               weight_matrix{i}(xpos,zpos) = w_norm;          
           end

           % plot weight matrix
           ax = subplot(rowNum,colNum, i);
           hold on;
           imagesc( [x_min x_max], [z_min z_max], weight_matrix{i}); 
           axis image;
           title(sprintf(['spotweights in beam ' num2str(i)]));
           xlabel('x [mm]');
           ylabel('z [mm]');
           cmap = colormap(ax, 'jet');
           cb = colorbar;
           ylabel(cb, 'normalized weight');
           hold off;
        end 
       
    elseif ( strcmp(stf(1).radiationMode, 'protons') || strcmp(stf(1).radiationMode, 'carbon') )
        
        % for protons and carbon Ions consider different energies as well
        % maximum weight
        w_max = max(weights);
        
        for i=1:numOfBeams
           fprintf(['Calculate weights for Beam ' num2str(i) '...\n']);
           
           % find range of ray positions within beam
           rayPos_mat = vertcat(stf(i).ray(:).rayPos_bev);
           x_min = min(rayPos_mat(:, 1));
           z_min = min(rayPos_mat(:, 3));
           x_max = max(rayPos_mat(:, 1));
           z_max = max(rayPos_mat(:, 3));

           % find all energies used
           all_energies{i} = unique(cat(2, stf(i).ray(:).energy));
           numOfEnergies(i) = length(all_energies{i});
           
           % initialise weight matrix
           weight_matrix{i} = zeros((x_max-x_min)/stf(i).bixelWidth,(z_max-z_min)/stf(i).bixelWidth, numOfEnergies(i));
           
           % find all bixel in this beam
           beam_ix = (bixelLut.beamNum == i);
           
           for j=1:numOfEnergies(i)
               % find all bixel with this energy
               energy_ix = and((bixelLut.energy == all_energies{i}(j)),beam_ix) ;
               
               for k=1:stf(i).numOfRays
                   % find weights for this ray
                   weight_ix = and((bixelLut.rayNum==k),energy_ix);

                   % normalize weight
                   w_norm = weights(weight_ix)/w_max;
                   if isempty(w_norm)
                       w_norm = 0;
                   end
                   
                   % find position in matrix corresponding to ray
                   xpos = (stf(i).ray(k).rayPos_bev(1) + abs(x_min))/stf(i).bixelWidth +1;
                   zpos = (stf(i).ray(k).rayPos_bev(3) + abs(z_min))/stf(i).bixelWidth +1;
                   weight_matrix{i}(xpos,zpos, j) = w_norm;          
               end
       
           end
           
           % new figure for each beam
           figure('Name', sprintf(['Beam ' num2str(i)]), 'WindowScrollWheelFcn', @figScroll);
           ax = axes;
           
           curr_ix_energy(i) = 1; % current energy slice index
           
           % plot spot weights
           img = imagesc( [x_min x_max], [z_min z_max], weight_matrix{i}(:,:,curr_ix_energy(i)));
           axis image;
           title(sprintf(['spotweights in beam ' num2str(i) ', energy: ' num2str(all_energies{i}(curr_ix_energy(i)))] ));
           xlabel('x [mm]');
           ylabel('z [mm]');
           cmap = colormap(ax, 'jet');
           caxis(ax, [0 1]);
           colormap(ax, cmap);
           cb = colorbar;
           ylabel(cb, 'normalized weight');
           hold off;  

        end

    end


    function figScroll(src,callbackdata)
      % function to enable scrolling through energy slices  
        
      % get handles for current figure and object
      curr_fig = gcf;
      curr_ax = gca(curr_fig);
      childrenOfAx = get(curr_ax,'Children');
      
      % get current beam number
      fig_name = get(curr_fig,'Name');
      curr_beam = str2double(fig_name(end));
      
      % set new energy slice
      curr_ix_energy(curr_beam) = curr_ix_energy(curr_beam) - callbackdata.VerticalScrollCount;
      
      % project to allowed range
      curr_ix_energy(curr_beam) = min(curr_ix_energy(curr_beam), numOfEnergies(curr_beam));
      curr_ix_energy(curr_beam) = max(curr_ix_energy(curr_beam), 1);
      
      % change Data in plot
      set(childrenOfAx,'CData',weight_matrix{curr_beam}(:,:,curr_ix_energy(curr_beam)));
      title(sprintf(['spotweights in beam ' num2str(curr_beam) ', energy: ' num2str(all_energies{curr_beam}(curr_ix_energy(curr_beam)))]));
      drawnow;
      
    end
end
