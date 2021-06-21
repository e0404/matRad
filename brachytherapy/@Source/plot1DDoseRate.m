function plot1DDoseRate(obj,r)
%PLOT1DDOSERATE Summary of this function goes here
%   Detailed explanation goes here
global cGy h mm

doseRate = obj.getDoseRate1D(r);

plot(r/mm,doseRate/(cGy/h))
xlabel('pos [mm]')
ylabel('Dose Rate cGy/h')


end

