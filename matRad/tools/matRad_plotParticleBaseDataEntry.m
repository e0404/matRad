function matRad_plotParticleBaseDataEntry(machine,index,hFigure)
%MATRAD_PLOTPARTICLEBASEDATAENTRY Summary of this function goes here
%   Detailed explanation goes here

if (ischar(machine) || isstring(machine)) && isfile(machine)
    load(machine);
end

if nargin < 3
    hFigure = figure;
end

subplot(2,3,1);
plot(machine.data(index).depths,machine.data(index).Z);
xlabel('depth [mm]');
ylabel('Z [MeV cm^2 /(g * primary)]');
title(sprintf('Depth Dose for E = %g MeV',machine.data(index).energy));

if isfield(machine.data(index),'sigma')
    subplot(2,3,2);
    plot(machine.data(index).depths,machine.data(index).sigma);
    xlabel('depth [mm]');
    ylabel('\sigma [mm]');
    title(sprintf('\\sigma for E = %g MeV',machine.data(index).energy));
end

if isfield(machine.data(index),'sigma')
    subplot(2,3,2);
    plot(machine.data(index).depths,machine.data(index).sigma); 
    xlabel('depth [mm]');
    ylabel('\sigma [mm]');
    title(sprintf('\\sigma for E = %g MeV',machine.data(index).energy));
end

if isfield(machine.data(index),'sigma1')
    subplot(2,3,2);
    plot(machine.data(index).depths,machine.data(index).sigma1);
    xlabel('depth [mm]');
    ylabel('\sigma [mm]');
    title(sprintf('\\sigma_1 for E = %g MeV',machine.data(index).energy));
end

if isfield(machine.data(index),'sigma2')
    subplot(2,3,3);
    plot(machine.data(index).depths,machine.data(index).sigma2);
    xlabel('depth [mm]');
    ylabel('\sigma [mm]');
    title(sprintf('\\sigma_2 for E = %g MeV',machine.data(index).energy));
end

if isfield(machine.data(index),'LET')
    subplot(2,3,4);
    plot(machine.data(index).depths,machine.data(index).LET);
    xlabel('depth [mm]');
    ylabel('LET [keV/\mu m]');
    title(sprintf('LET for E = %g MeV',machine.data(index).energy));
end

if isfield(machine.data(index),'alpha')
    subplot(2,3,5);
    nT = numel(machine.data(index).alphaBetaRatio);
    legendNames = cell(1,nT);
    for i =1:nT
        plot(machine.data(index).depths,machine.data(index).alpha); hold on;
        legendNames{i} = sprintf('[\\alpha_\\gamma = %g, \\beta_\\gamma = %g]',machine.data(index).alphaX(i),machine.data(index).betaX(i));
    end
    legend(legendNames);

    xlabel('depth [mm]');
    ylabel('\alpha [Gy^{-1}]');
    title(sprintf('\\alpha for E = %g MeV',machine.data(index).energy));
end

if isfield(machine.data(index),'beta')
    subplot(2,3,6);
    nT = numel(machine.data(index).alphaBetaRatio);
    legendNames = cell(1,nT);
    for i =1:nT
        plot(machine.data(index).depths,machine.data(index).beta); hold on;
        legendNames{i} = sprintf('[\\alpha_\\gamma = %g, \\beta_\\gamma = %g]',machine.data(index).alphaX(i),machine.data(index).betaX(i));
    end
    legend(legendNames);

    xlabel('depth [mm]');
    ylabel('\beta [Gy^{-2}]');
    title(sprintf('\\beta for E = %g MeV',machine.data(index).energy));
end

end

