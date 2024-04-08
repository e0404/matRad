function [density, sound_speed] = hounsfield2density(ct_data, plot_fitting)
%HOUNSFIELD2DENSITY Convert Hounsfield units to density.
%
% DESCRIPTION:
%     hounsfield2density converts Hounsfield units to units of density 
%     [kg/m^3] based on the experimental data given by Schneider et al. The
%     conversion is made using a piece-wise linear fit to the data. For
%     soft-tissue, the approximate sound speed can also be returned using
%     the empirical relationship given by Mast. 
%
% USAGE:
%     density = hounsfield2density(ct_data)
%     density = hounsfield2density(ct_data, plot_fitting)
%     [density, sound_speed] = hounsfield2density(ct_data)
%     [density, sound_speed] = hounsfield2density(ct_data, plot_fitting)
%
% INPUTS:
%     ct_data      - CT data in Hounsfield units to convert to density
%
% OPTIONAL INPUTS:
%     plot_fitting - Boolean controlling whether the original data points
%                    and fitting is plotted (default = false)
%
% OUTPUTS:
%     density      - density in kg/m^3
%     sound_speed  - sound speed in m/s
%
% ABOUT:
%     author       - Bradley Treeby
%     date         - 9th January 2012
%     last update  - 4th June 2017
%
% REFERENCES:
%     Schneider, U., Pedroni, E., and Lomax A., "The calibration of CT
%     Hounsfield units for radiotherapy treatment planning," Phys. Med.
%     Biol., 41, pp. 111-124 (1996).
%
%     Mast, T. D., "Empirical relationships between acoustic parameters
%     in human soft tissues," Acoust. Res. Lett. Online, 1(2), pp. 37-42
%     (2000). 
%
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2012-2017 Bradley Treeby

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>. 

% create empty density matrix
density = zeros(size(ct_data), 'like', ct_data);

% apply conversion in several parts using linear fits to the data
% Part 1: Less than 930 Hounsfield Units
density(ct_data < 930) = polyval([1.025793065681423, -5.680404011488714], ct_data(ct_data < 930));

% Part 2: Between 930 and 1098 (soft tissue region)
density(ct_data >= 930 & ct_data <= 1098) = polyval([0.9082709691264, 103.6151457847139], ct_data(ct_data >= 930 & ct_data <= 1098));

% Part 3: Between 1098 and 1260 (between soft tissue and bone)
density(ct_data > 1098 & ct_data < 1260) = polyval([0.5108369316599, 539.9977189228704], ct_data(ct_data > 1098 & ct_data < 1260));

% Part 4: Greater than 1260 (bone region)
density(ct_data >= 1260) = polyval([0.6625370912451, 348.8555178455294], ct_data(ct_data >= 1260));

% calculate corresponding sound speed values if required using soft tissue
% relationship
if nargout == 2
    sound_speed = (density + 349) ./ 0.893;
end

% plot original data and fitted curves if required
if nargin == 2 && plot_fitting
        
    % soft tissue values excluding spongiosa
    density_soft_tissue = [0.95, 1.06, 1.04, 1.02, 1.00, 1.07, 1.03, 1.06, 1.05, 1.06, 1.05, 1.03, 1.05, 1.05, 1.04, 1.10, 1.03, 0.98, 1.09, 1.06, 1.04, 1.05] * 1000;
    hounsfd_soft_tissue = [ 930, 1055, 1037, 1003, 1003, 1050, 1023, 1055, 1043, 1053, 1044, 1028, 1042, 1045, 1032, 1098, 1014,  958, 1075, 1054, 1032, 1040];

    % bone values
    density_bone = [1.92, 1.61, 1.33, 1.46, 1.68, 1.41, 1.52, 1.29, 1.18, 1.42, 1.33] * 1000;
    hounsfd_bone = [2376, 1903, 1499, 1683, 2006, 1595, 1763, 1413, 1260, 1609, 1477];
    
    % filled lung values
    density_lung = 0.26 * 1000;
    hounsfd_lung = 259;

    % find linear fit for soft tissue data points
    h_axis_soft_tissue = min(hounsfd_soft_tissue(:)):max(hounsfd_soft_tissue(:));
    p_soft_tissue = polyfit(hounsfd_soft_tissue, density_soft_tissue, 1);
    
    % find linear fit for bone data points
    h_axis_bone = min(hounsfd_bone(:)):max(hounsfd_bone(:));
    p_bone = polyfit(hounsfd_bone, density_bone, 1);

    % find linear fit from soft tissue to bone region
    h_axis_tissue_to_bone = [h_axis_soft_tissue(end), h_axis_bone(1)];
    density_tissue_to_bone = [polyval(p_soft_tissue, h_axis_soft_tissue(end)), polyval(p_bone, h_axis_bone(1))];
    p_tissue_to_bone = polyfit(h_axis_tissue_to_bone, density_tissue_to_bone, 1);
    
    % find linear fit from filled lung to soft tissue region
    h_axis_lung_to_tissue = [hounsfd_lung, h_axis_soft_tissue(1)];
    density_lung_to_tissue = [density_lung, polyval(p_soft_tissue, h_axis_soft_tissue(1))];
    p_lung_to_tissue = polyfit(h_axis_lung_to_tissue, density_lung_to_tissue, 1);
    
    % plot original data points
    figure;
    hold on;
    plot(hounsfd_soft_tissue, density_soft_tissue, 'r.');
    plot(hounsfd_bone , density_bone , 'b.');
    plot(hounsfd_lung, density_lung, 'g.');
    
    % plot fitting data
    plot(h_axis_soft_tissue, polyval(p_soft_tissue, h_axis_soft_tissue), 'r-');
    plot(h_axis_bone, polyval(p_bone, h_axis_bone), 'b-');
    plot(h_axis_tissue_to_bone, polyval(p_tissue_to_bone, h_axis_tissue_to_bone), 'k-');
    plot(h_axis_lung_to_tissue, polyval(p_lung_to_tissue, h_axis_lung_to_tissue), 'k-');
    xlabel('Hounsfield Units');
    ylabel('Density [kg/m^3]');
    legend('Soft Tissue', 'Bone', 'Lung', 'Location', 'NorthWest');
    box on;

end