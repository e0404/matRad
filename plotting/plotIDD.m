function plotIDD(varargin)
%PLOTIDD Summary of this function goes here
%   Detailed explanation goes here
figure
hold on

if ~isstruct(varargin{end-1})
    profile = varargin{end-1};
    until = numel(varargin)-2;
else
    profile = 0;
    until = numel(varargin)-1;
end

if isstr(varargin{end}) && (strcmp(varargin{end},'RBE') || strcmp(varargin{end},'RBExD'))
    for i = 1:until
        plot(matRad_calcIDD(varargin{i}.RBExD,profile),'LineWidth',1.5, 'DisplayName', inputname(i))
    end
if profile
    ylabel('profile / RBExD')
 else   
    ylabel('IDD / RBExD')
    end
else
    for i = 1:until+1
        plot(matRad_calcIDD(varargin{i}.physicalDose,profile),'LineWidth',1.5, 'DisplayName', inputname(i))
    end
    if profile
    ylabel('profile / physicalDose')
 else   
    ylabel('IDD / physicalDose')
    end
end
xlabel('depth [voxel]')

% xlim([0 size(varargin{i}.physicalDose,1)*0.8])
legend('Interpreter', 'none','Location','northwest')
end

