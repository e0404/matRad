% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Decription: Generate time array formatted for log file creation
%
% Author: Lucas Sommer (Lucas.Sommer@tum.de), 10/2018
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function timeArray = matRad_getTime4log
uhrzeit = datestr(now,'HH:MM:SS');
datum  =  datestr(now,'yyyy-mm-dd');

timeArray = strcat(datum(1:4), datum(6:7), datum(9:10), 'd', uhrzeit(1:2), ...
    uhrzeit(4:5), uhrzeit(7:8), 'h');