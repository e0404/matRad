function bixelInfo = matRad_makeBixelTimeSeq(stf, resultGUI, plotting)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
% call
%   
%
% input
%       
%  
% output
%
% comment:
% 
% References
%
% Ahmad Neishabouri June 2018
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    plotting = 'off';
end

% defining the constant parameters
%
% time required for synchrotron to change energy
es_time = 4 * 10^6;
% time required for synchrotron to recharge it's spill
spill_recharge_time = 2 * 10^6;
% number of particles generated in each spill
spill_size = 4 * 10 ^ 10;
% speed of synchrotron's lateral scanning in an IES
scan_speed = 10; % m/s
% number of particles per second
spill_intensity = 4 * 10 ^ 8;


for i = 1:length(stf)
    steerTime(i) = stf(i).bixelWidth * (10 ^ 3)/ scan_speed;
end

bixelInfo = struct;

% first loop loops over all bixels to store their position and ray number
% in each IES
wOffset = 0;
for i = 1:length(stf) % looping over all beams
    
    usedEnergies = unique([stf(i).ray(:).energy]);
    usedEnergiesSorted = sort(usedEnergies, 'descend');
    
    bixelInfo(i).orderToSTF = zeros(stf(i).totalNumOfBixels, 1);
    bixelInfo(i).time = zeros(stf(i).totalNumOfBixels, 1);
    bixelInfo(i).e = zeros(stf(i).totalNumOfBixels, 1);
    bixelInfo(i).orderToSS = zeros(stf(i).totalNumOfBixels, 1);
    
    for e = 1:length(usedEnergies) % looping over IES's
        
        s = 1;
        
        for j = 1:stf(i).numOfRays % looping over all rays
            
            % find the rays which are active in current IES
            if(any(stf(i).ray(j).energy == usedEnergiesSorted(e)))
                
                x = stf(i).ray(j).rayPos_bev(1);
                y = stf(i).ray(j).rayPos_bev(3);
                %
                bixelInfo(i).IES(e).x(s)       = x; % store x position
                bixelInfo(i).IES(e).y(s)       = y; % store y position
                bixelInfo(i).IES(e).w_index(s) = wOffset + ...
                    sum(stf(i).numOfBixelsPerRay(1:(j-1))) + ...
                    find(stf(i).ray(j).energy == usedEnergiesSorted(e)); % store index
                
                s = s + 1;
                
            end
        end
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
    
    temp = 0; % temporary variable for plotting figures
    t = 0;
    
    for e = 1: length(usedEnergies)
        
        % sort the y positions from high to low (backforth is up do down)
        y_sorted = sort(unique(bixelInfo(i).IES(e).y), 'descend');
        x_sorted = sort(bixelInfo(i).IES(e).x, 'ascend');
        
        for k = 1:length(y_sorted)
            
            y = y_sorted(k);
            % find indexes corresponding to current y position
            % in other words, number of bixels in the current row
            ind_y = find(bixelInfo(i).IES(e).y == y);
            
            % since backforth fasion is zig zag like, flip the order every
            % second row
            if ~rem(k,2)
                ind_y = fliplr(ind_y);
            end
            
            % loop over all the bixels in the row
            for s = ind_y
                
                x = bixelInfo(i).IES(e).x(s);
                % TODO: sort x's in case someone hacked stf
                
                w_index = bixelInfo(i).IES(e).w_index(s);
                
                % calculating the time:
                %
                % required spot fluence
                protons = resultGUI.w(w_index)* 10^6;
                % time spent to spill the required spot fluence
                spillTime = protons * 10^6 / spill_intensity;
                %
                % spotTime:time spent to steer scan along IES per bixel
                t = t + steerTime(i) + spillTime;
                %
                % taking account of the time to recharge the spill in case
                % the required fluence was more than spill size
                if(spill_usage + protons > spill_size)
                    t = t + spill_recharge_time;
                    spill_usage = 0;
                end
                %
                % used amount of fluence from current spill
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
                
                % storing the time and the order of bixels
                %
                % make the both counter and index 'per beam' - help index
                order_count = order_counter - offset;
                w_ind = w_index - offset;
                %
                % timeline according to the spot scanning order
                bixelInfo(i).time(order_count) = t;
                % IES of bixels according to the spot scanning order
                bixelInfo(i).e(order_count) = e;
                % according to spot scanning order, sorts w index of all
                % bixels, use this order to transfer STF order to Spot
                % Scanning order
                bixelInfo(i).orderToSS(order_count) = w_ind;
                %
                % according to STF order, gives us order of irradiation of
                % each bixel, use this order to transfer Spot Scanning
                % order to STF order
                % orderToSTF(orderToSS) = orderToSS(orderToSTF) = 1:#bixels
                bixelInfo(i).orderToSTF(w_ind) = order_count;
                %
                order_counter  = order_counter + 1;
                
            end
        end
        
        t = t + es_time;
        
    end
    
    % storing the fluence per beam
    bixelInfo(i).w = resultGUI.w(offset + 1: offset + stf(i).totalNumOfBixels);
    
    offset = offset + stf(i).totalNumOfBixels;
    
end

end