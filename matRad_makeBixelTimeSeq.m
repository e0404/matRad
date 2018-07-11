function bixelInfo = matRad_makeBixelTimeSeq(stf, w, plotting)


if nargin < 3
    plotting = 'off';
end

mil = 10^6;
es_time = 4 * mil;
spill_recharge_time = 2 * mil;
scan_speed = 10; % m/s
spill_size = 4 * 10 ^ 10;
spill_intensity = 4 * 10 ^ 8; 

for i = 1:length(stf)
    spot_time(i) = stf(i).bixelWidth * (10 ^ 3)/ scan_speed;
end

order = zeros(sum([stf.totalNumOfBixels]), 1);
bixelInfo = struct;

% first loop loops over all bixels to store their position and ray number
% in each IES
wOffset = 0;
for i = 1:length(stf) % looping over all beams
    
    usedEnergies = unique([stf(i).ray(:).energy]);
    usedEnergiesSorted = sort(usedEnergies, 'descend');
    
    bixelInfo(i).order = zeros(stf(i).totalNumOfBixels, 1);
    bixelInfo(i).time = zeros(stf(i).totalNumOfBixels, 1);
    bixelInfo(i).e = zeros(stf(i).totalNumOfBixels, 1);
    
    bixelInfo(i).w_index = [];
    
    for e = 1:length(usedEnergies) % looping over IES's
        
        s = 1;
        
        for j = 1:stf(i).numOfRays % looping over all rays
            
            % find the rays which are active in current IES
            if(any(stf(i).ray(j).energy == usedEnergiesSorted(e)))
                
                x = stf(i).ray(j).rayPos_bev(1);
                y = stf(i).ray(j).rayPos_bev(3);

                bixelInfo(i).IES(e).x(s)       = x; % store x position
                bixelInfo(i).IES(e).y(s)       = y; % store y position
                bixelInfo(i).IES(e).w_index(s) = wOffset + ...
                                                 sum(stf(i).numOfBixelsPerRay(1:(j-1))) + ...
                                                 find(stf(i).ray(j).energy == usedEnergiesSorted(e)); % store index
                
                s = s + 1;
                
            end
        end
        
    bixelInfo(i).w_index = [bixelInfo(i).w_index; bixelInfo(i).IES(e).w_index'];
        
    end
    
    wOffset = wOffset + sum(stf(i).numOfBixelsPerRay);
    
end

%
% after storing all the required information,
% same loop over all bixels will put each bixel in it's order

order_counter = 1;
spill_usage = 0;
offset = 0;

for i = 1:length(stf)
    
    usedEnergies = unique([stf(i).ray(:).energy]);
    usedEnergiesSorted = sort(usedEnergies, 'descend');
    temp = 0; % temporary variable for plotting figures
    t = 0;
    
    for e = 1: length(usedEnergies)
        
        % sort the y positions from high to low (backforth is up do down)
        y_sorted = sort(unique(bixelInfo(i).IES(e).y), 'descend');
        x_sorted = sort(bixelInfo(i).IES(e).x, 'ascend');
        
        for k = 1:length(y_sorted)
            
            y = y_sorted(k);
            x = x_sorted(k);
            % find indexes corresponding to current y position
            % in other words, number of bixels in the current row
            index = find(bixelInfo(i).IES(e).y == y);
            
            % since backforth fasion is zig zag like, flip the order every
            % second row
            if ~rem(k,2)
                index = fliplr(index);
            end
            
            % loop over all the bixels in the row
            for ss = index
                
                x = bixelInfo(i).IES(e).x(ss);
                
                w_index = bixelInfo(i).IES(e).w_index(ss);
                
                
                protons = w(w_index)* 10^6;
                
                spill_time = protons * 10^6 / spill_intensity;
                
                t = t + spot_time(i) + spill_time;
                
                if(spill_usage + protons > spill_size)
                    t = t + spill_recharge_time;
                    spill_usage = 0;
                end
                
                spill_usage = spill_usage + protons;
                
                % following is to plot the bixel ordering, can be deleted
                % in the final version of the code
                if(strcmp(plotting, 'on'))
                    
                    if(temp ~= e)
                        clf
                        h = animatedline('LineStyle', 'none', 'Marker', 'o');
                        axis([-50 50 -50 50])
                        title(['Beam #', num2str(i), ', IES #', num2str(e)])
                    end
                    
                    temp = e;
                    
                    addpoints(h, x, y);
                    strmin = [num2str(order_counter), '  '];
                    text(x,y,strmin,'HorizontalAlignment','right');
                    pause(.1);
                    drawnow
                end
                
                % temp solution!
                order_ind = order_counter - offset;
                w_ind = w_index - offset;
                
                % assign the order to the corresponding stf index
                bixelInfo(i).time(order_ind) = t;
                bixelInfo(i).order(w_ind) = order_ind;
                bixelInfo(i).w(w_ind) = w(w_index);
                bixelInfo(i).e(order_ind) = e;
                
                order(w_index) = order_counter;
                order_counter  = order_counter + 1;
            end
        end
        
        t = t + es_time;
    end
    
    offset = offset + stf(i).totalNumOfBixels;
    
    
end

end
%}