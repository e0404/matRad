function MCsettings = matRad_MCinit(pln,MCsettings)

if ~isfield(MCsettings,'MC_PBS')
    MCsettings.MC_PBS = 0;
end

if ~isfield(MCsettings,'MCengine')
    if strcmp(pln.radiationMode,'protons')
        MCsettings.MCengine = 'MCsquare';
    else
        error('No MC engine available for other particles yet!');
    end
end

if ~isfield(MCsettings,'runsPath')
    MCsettings.runsPath = 'MCexport';
end

if ~isfield(MCsettings,'fractionHistories')
    MCsettings.fractionHistories = 1e-4;
end

if ~isfield(MCsettings,'ElectronProdCut_mm')
    MCsettings.runsPath = 0.5;
end

if ~isfield(MCsettings,'minRelWeight')
    MCsettings.minRelWeight = 0;
end

end