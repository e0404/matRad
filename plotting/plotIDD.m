function plotIDD(varargin)
%PLOTIDD Summary of this function goes here
%   Detailed explanation goes here
figure
hold on

if isstr(varargin{end}) && (strcmp(varargin{end},'RBE') || strcmp(varargin{end},'RBExD'))
    for i = 1:numel(varargin)-1
        plot(matRad_calcIDD(varargin{i}.RBExD,1),'LineWidth',1.5, 'DisplayName', inputname(i))
    end
    ylabel('RBExD')
else
    for i = 1:numel(varargin)
        plot(matRad_calcIDD(varargin{i}.physicalDose,1),'LineWidth',1.5, 'DisplayName', inputname(i))
    end
    ylabel('physicalDose')
end
xlabel('depth [voxel]')
% xlim([0 size(varargin{i}.physicalDose,1)*0.8])
legend('Interpreter', 'none','Location','northwest')
end

