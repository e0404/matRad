function timeSequence = matRad_makeBixelTimeSeq(stf, resultGUI)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% using the steering information of matRad, makes a time sequenced order
% according to the irradiation scheme in spot scanning
%
% call
%   timeSequence = matRad_makeBixelTimeSeq(stf, resultGUI)
%
% input
%   stf:            matRad steering information struct
%   resultGUI:      struct containing optimized fluence vector
%
% output
%   timeSequence:      struct containing bixel ordering information and the
%                   time sequence of the spot scanning
%
% References
%   spill structure and timing informations:
%   http://cdsweb.cern.ch/record/1182954
%   http://iopscience.iop.org/article/10.1088/0031-9155/56/20/003/meta
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2018 the matRad development team.
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% defining the constant parameters
%
% time required for synchrotron to change energy

es_time = 3 * 10^6; % [\mu s]
% time required for synchrotron to recharge it' spill
spill_recharge_time = 2 * 10^6; % [\mu s]
% number of particles generated in each spill
spill_size = 4 * 10 ^ 10;
% speed of synchrotron's lateral scanning in an IES
scan_speed = 10; % m/s
% number of particles per second
spill_intensity = 4 * 10 ^ 8;


steerTime = [stf.bixelWidth] * (10 ^ 3)/ scan_speed; % [\mu s]

timeSequence = struct;

% first loop loops over all bixels to store their position and ray number
% in each IES
wOffset = 0;
for i = 1:length(stf) % looping over all beams
    
    usedEnergies = unique([stf(i).ray(:).energy]);
    usedEnergiesSorted = sort(usedEnergies, 'descend');
    
    timeSequence(i).orderToSTF = zeros(stf(i).totalNumOfBixels, 1);
    timeSequence(i).orderToSS = zeros(stf(i).totalNumOfBixels, 1);
    timeSequence(i).time = zeros(stf(i).totalNumOfBixels, 1);
    timeSequence(i).e = zeros(stf(i).totalNumOfBixels, 1);
    
    
    for e = 1:length(usedEnergies) % looping over IES's
        
        s = 1;
        
        for j = 1:stf(i).numOfRays % looping over all rays
            
            % find the rays which are active in current IES
            if(any(stf(i).ray(j).energy == usedEnergiesSorted(e)))
                
                x = stf(i).ray(j).rayPos_bev(1);
                y = stf(i).ray(j).rayPos_bev(3);
                %
                timeSequence(i).IES(e).x(s)       = x; % store x position
                timeSequence(i).IES(e).y(s)       = y; % store y position
                timeSequence(i).IES(e).w_index(s) = wOffset + ...
                    sum(stf(i).numOfBixelsPerRay(1:(j-1))) + ...
                    find(stf(i).ray(j).energy == usedEnergiesSorted(e)); % store index
                
                s = s + 1;
                
            end
        end
    end
    
    wOffset = wOffset + sum(stf(i).numOfBixelsPerRay);
    
end

% after storing all the required information,
% same loop over all bixels will put each bixel in it's order

spill_usage = 0;
offset = 0;

for i = 1:length(stf)
    
    usedEnergies = unique([stf(i).ray(:).energy]);
    
    t = 0;
    order_count = 1;
    
    for e = 1: length(usedEnergies)
        
        % sort the y positions from high to low (backforth is up do down)
        y_sorted = sort(unique(timeSequence(i).IES(e).y), 'descend');
        x_sorted = sort(timeSequence(i).IES(e).x, 'ascend');
        
        for k = 1:length(y_sorted)
            
            y = y_sorted(k);
            % find indexes corresponding to current y position
            % in other words, number of bixels in the current row
            ind_y = find(timeSequence(i).IES(e).y == y);
            
            % since backforth fasion is zig zag like, flip the order every
            % second row
            if ~rem(k,2)
                ind_y = fliplr(ind_y);
            end
                        
            % loop over all the bixels in the row
            for is = 1:length(ind_y)
                
                s = ind_y(is);
                
                x = x_sorted(s);
                
                w_index = timeSequence(i).IES(e).w_index(s);
                
                % in case there were holes inside the plan "multi"
                % multiplies the steertime to take it into account:
                if(k == 1 && is == 1)
                    x_prev = x;
                    y_prev = y;
                end
                % x direction
                multi = abs(x_prev - x)/stf(i).bixelWidth;
                % y direction
                multi = multi + abs(y_prev - y)/stf(i).bixelWidth;
                %
                x_prev = x;
                y_prev = y;
                
                % calculating the time:
                
                % required spot fluence
                numOfParticles = resultGUI.w(w_index)* 10^6;
                % time spent to spill the required spot fluence
                spillTime = numOfParticles * 10^6 / spill_intensity;
                
                % spotTime:time spent to steer scan along IES per bixel
                t = t + multi * steerTime(i) + spillTime;

                % taking account of the time to recharge the spill in case
                % the required fluence was more than spill size
                if(spill_usage + numOfParticles > spill_size)
                    t = t + spill_recharge_time;
                    spill_usage = 0;
                end
                
                % used amount of fluence from current spill
                spill_usage = spill_usage + numOfParticles;
                
                % storing the time and the order of bixels
                
                % make the both counter and index 'per beam' - help index
                w_ind = w_index - offset;
                
                % timeline according to the spot scanning order
                timeSequence(i).time(order_count) = t;
                % IES of bixels according to the spot scanning order
                timeSequence(i).e(order_count) = e;
                % according to spot scanning order, sorts w index of all
                % bixels, use this order to transfer STF order to Spot
                % Scanning order
                timeSequence(i).orderToSS(order_count) = w_ind;

                % according to STF order, gives us order of irradiation of
                % each bixel, use this order to transfer Spot Scanning
                % order to STF order
                % orderToSTF(orderToSS) = orderToSS(orderToSTF) = 1:#bixels
                timeSequence(i).orderToSTF(w_ind) = order_count;
                
                order_count  = order_count + 1;
                
            end
        end
        
        t = t + es_time;
        
    end
    
    % storing the fluence per beam
    timeSequence(i).w = resultGUI.w(offset + 1: offset + stf(i).totalNumOfBixels);
    
    offset = offset + stf(i).totalNumOfBixels;
    
end

end