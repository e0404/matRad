function extractPlans(obj)
%EXTRACTPLANS Summary of this function goes here
%   Detailed explanation goes here
global seconds
for ip=1:length(obj.RPfiles)
    inf     = dicominfo(obj.RPfiles{ip});
    actItem = ['Item_',num2str(ip)];
    
    if isfield(inf,'ApplicationSetupSequence')
        config  = fieldnames(inf.ApplicationSetupSequence);
        
        if length(config)>1
            %%
            noSeeds = length(config);
            obj.RP.(actItem).seeds.dwellPosition = [];
            for is =1:noSeeds
                
                if isfield(inf.ApplicationSetupSequence.(config{is}),'ChannelSequence')
                    cNames = fieldnames(inf.ApplicationSetupSequence.(config{is}).ChannelSequence);
                    
                    if length(cNames)==1
                        if inf.ApplicationSetupSequence.(config{is}).ChannelSequence.(cNames{1}).BrachyControlPointSequence.Item_2.CumulativeTimeWeight==100
                            obj.RP.(actItem).seeds.dwellPosition(end+1,:) = inf.ApplicationSetupSequence.(config{is}).ChannelSequence.(cNames{1}).BrachyControlPointSequence.Item_2.ControlPoint3DPosition;
                        else
                            error('extractPlan.m: Cumulative Time Weight should be 100')
                        end
                    else
                        disp('extractPlan.m: only one item for ChannelSequence is expected')
                    end
                end
            end
            disp('extractPlan.m: more than one fields in ApplicationSetupSequence found (i.e. plan without needles)')
        else
            for ic = 1:length(config)
                
                needles = fieldnames(inf.ApplicationSetupSequence.(config{ic}).ChannelSequence);
                obj.RP.(actItem).additionalInfo.finalCumTimeWeight = zeros(length(needles),1);
                obj.RP.(actItem).additionalInfo.channelTotTime     = zeros(length(needles),1);
                obj.RP.(actItem).additionalInfo.channelLength      = zeros(length(needles),1);
                obj.RP.(actItem).additionalInfo.stepSize           = zeros(length(needles),1);
                obj.RP.(actItem).additionalInfo.templateRowCol     = cell(length(needles),1);
                if isfield(inf.ApplicationSetupSequence.Item_1,'Private_300b_1011')
                    obj.RP.(actItem).additionalInfo.tandemInfo  = inf.ApplicationSetupSequence.Item_1.Private_300b_1011;
                end
                obj.RP.(actItem).seeds.dwellPosition = cell(length(needles),1);
                obj.RP.(actItem).seeds.dwellTime     = cell(length(needles),1);
                obj.RP.(actItem).seeds.dwellWeight   = cell(length(needles),1);
                obj.RP.(actItem).seeds.name          = cell(length(needles),1);
                
                for in = 1:length(needles)
                    if isfield(inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}),'SourceApplicatorType')
                        obj.RP.(actItem).seeds.name{in}                        = inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}).SourceApplicatorType;
                    end
                    if isfield(inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}),'ChannelLength')
                        obj.RP.(actItem).additionalInfo.channelLength(in)      = inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}).ChannelLength;
                    end
                    if isfield(inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}),'SourceApplicatorStepSize')
                        obj.RP.(actItem).additionalInfo.stepSize(in)           = inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}).SourceApplicatorStepSize;
                    end
                    obj.RP.(actItem).additionalInfo.channelTotTime(in)     = inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}).ChannelTotalTime;
                    if isfield(inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}),'FinalCumulativeTimeWeight')
                        obj.RP.(actItem).additionalInfo.finalCumTimeWeight(in) = inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}).FinalCumulativeTimeWeight;
                    end
                    
                    if isfield(inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}),'Private_1001_10xx_Creator')
                        obj.RP.(actItem).additionalInfo.templateRowCol{in} = double([inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}).Private_1001_107a,inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}).Private_1001_107b]);
                    end
                    
                    if isfield(inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}),'Private_1001_1080')
                        nItem = fieldnames(inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}).Private_1001_1080);
                        pos3D = zeros(length(nItem),3);
                        for ini=1:length(nItem)
                            pos3D(ini,:) = inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}).Private_1001_1080.(nItem{ini}).Private_1001_1083';
                        end
                        
                        obj.RP.(actItem).needles.pos{in} = pos3D;
                    end
                       
                    
                    if isfield(inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}),'BrachyControlPointSequence');
                        posItem = fieldnames(inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}).BrachyControlPointSequence);
                        pos3D   = zeros(length(posItem),3);
                        relPos  = zeros(length(posItem),1);
                        contP   = zeros(length(posItem),1);
                        cumTime = zeros(length(posItem),1);
                        
                        %%
                        pos3DFound = true;
                        for ipi=1:length(posItem)
                            if isfield(inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}).BrachyControlPointSequence.(posItem{ipi}),'ControlPoint3DPosition')
                                pos3D(ipi,:) = inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}).BrachyControlPointSequence.(posItem{ipi}).ControlPoint3DPosition';
                                
                                if isempty(inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}).BrachyControlPointSequence.(posItem{ipi}).CumulativeTimeWeight)
                                    cumTime(ipi) = 0;
                                else
                                    cumTime(ipi) = inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}).BrachyControlPointSequence.(posItem{ipi}).CumulativeTimeWeight;
                                end
                                relPos(ipi)  = inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}).BrachyControlPointSequence.(posItem{ipi}).ControlPointRelativePosition;
                                contP(ipi)   = inf.ApplicationSetupSequence.(config{ic}).ChannelSequence.(needles{in}).BrachyControlPointSequence.(posItem{ipi}).ControlPointIndex;
                            else
                                pos3DFound = false;
                            end
                        end
                        %%
                        obj.RP.(actItem).seeds.needleID(in)      = in;
                        if pos3DFound
                            
                            time = obj.RP.(actItem).additionalInfo.channelTotTime(in)*cumTime./obj.RP.(actItem).additionalInfo.finalCumTimeWeight(in);
                            tempTime = zeros(length(posItem)/2,1);
                            tempPos  = zeros(length(posItem)/2,3);
                            tempWeight = zeros(length(posItem)/2,1);
                            for i=2:2:length(posItem)
                                tempTime(i/2)   = time(i)-time(i-1);
                                tempPos(i/2,:)  = pos3D(i,:);
                                tempWeight(i/2) = cumTime(i)-cumTime(i-1);
                            end
                            

                            
                            obj.RP.(actItem).seeds.dwellPosition{in} = tempPos;
                            if obj.RP.(actItem).additionalInfo.finalCumTimeWeight(in)>0
                                obj.RP.(actItem).seeds.dwellWeight{in}   = tempWeight;
                                obj.RP.(actItem).seeds.dwellTime{in}     = tempTime*seconds;
                            else
                                obj.RP.(actItem).seeds.dwellWeight{in} = zeros(size(tempPos,1),1);
                                obj.RP.(actItem).seeds.dwellTime{in}   = zeros(size(tempPos,1),1);
                            end
                        end
                    end
                end
                
                indActive = ~cellfun(@isempty,obj.RP.(actItem).seeds.dwellPosition);
                obj.RP.(actItem).seeds.dwellPosition = obj.RP.(actItem).seeds.dwellPosition(indActive);
                obj.RP.(actItem).seeds.dwellTime     = obj.RP.(actItem).seeds.dwellTime(indActive);
                obj.RP.(actItem).seeds.dwellWeight   = obj.RP.(actItem).seeds.dwellWeight(indActive);
                obj.RP.(actItem).seeds.name          = obj.RP.(actItem).seeds.name(indActive);
                obj.RP.(actItem).seeds.needleID      = obj.RP.(actItem).seeds.needleID(indActive);
                
                obj.RP.(actItem).additionalInfo.finalCumTimeWeight = obj.RP.(actItem).additionalInfo.finalCumTimeWeight(indActive);
                obj.RP.(actItem).additionalInfo.channelTotTime     = obj.RP.(actItem).additionalInfo.channelTotTime(indActive);
                obj.RP.(actItem).additionalInfo.channelLength      = obj.RP.(actItem).additionalInfo.channelLength(indActive);
                obj.RP.(actItem).additionalInfo.stepSize           = obj.RP.(actItem).additionalInfo.stepSize(indActive);
                obj.RP.(actItem).additionalInfo.templateRowCol     = obj.RP.(actItem).additionalInfo.templateRowCol(indActive);
                
                if isfield(inf,'DoseReferenceSequence')
                    refPoints = fieldnames(inf.DoseReferenceSequence);
                    obj.RP.(actItem).refPoints.position = zeros(length(refPoints),3);
                    obj.RP.(actItem).refPoints.dose     = zeros(length(refPoints),1);
                    obj.RP.(actItem).refPoints.name     = cell(length(refPoints),1);
                    for irp=1:length(refPoints)
                        obj.RP.(actItem).refPoints.position(irp,:) = inf.DoseReferenceSequence.(refPoints{irp}).DoseReferencePointCoordinates';
                        obj.RP.(actItem).refPoints.dose(irp)       = inf.DoseReferenceSequence.(refPoints{irp}).TargetPrescriptionDose;
                        obj.RP.(actItem).refPoints.name{irp}       = inf.DoseReferenceSequence.(refPoints{irp}).DoseReferenceDescription;
                    end
                else
                    obj.RP.(actItem).refPoints = [];
                end
                
                
                tempTime   = cell2mat(obj.RP.(actItem).seeds.dwellTime);
                tempWeight = cell2mat(obj.RP.(actItem).seeds.dwellWeight);
                temp       = tempWeight./tempTime;
                if std(temp(~isnan(temp)))/mean(temp(~isnan(temp)))<1e-6
                    obj.RP.(actItem).additionalInfo.time2Weight = mean(temp(~isnan(temp)));
                else
                    obj.RP.(actItem).additionalInfo.time2Weight = [];
                end
            end
        end
        if isfield(inf,'StudyDescription')
            obj.RP.(actItem).additionalInfo.description     = inf.StudyDescription;
        else
            obj.RP.(actItem).additionalInfo.description     = [];
        end
        if isfield(inf,'ApprovalStatus')
            obj.RP.(actItem).additionalInfo.status          = inf.ApprovalStatus;
        else
            obj.RP.(actItem).additionalInfo.status          = [];
        end
        if isfield(inf,'SourceStrengthUnits')
            obj.RP.(actItem).source.SourceStrengthUnits     = inf.SourceSequence.Item_1.SourceStrengthUnits;
        else
            obj.RP.(actItem).source.SourceStrengthUnits     = [];
        end
        
        if isfield(inf,'ReferencedStructureSetSequence')
        obj.RP.(actItem).additionalInfo.structID            = inf.ReferencedStructureSetSequence.Item_1.ReferencedSOPInstanceUID;
        else
            obj.RP.(actItem).additionalInfo.structID            = [];
        end
        obj.RP.(actItem).additionalInfo.prescribedDose      = inf.FractionGroupSequence.Item_1.ReferencedBrachyApplicationSetupSequence.Item_1.BrachyApplicationSetupDose;
        obj.RP.(actItem).additionalInfo.BrachyTreatmentType = inf.BrachyTreatmentType;
        obj.RP.(actItem).additionalInfo.filename            = obj.RPfiles{ip};
        obj.RP.(actItem).source.Isotope                     = inf.SourceSequence.Item_1.SourceIsotopeName;
        obj.RP.(actItem).source.SourceIsotopeHalfLife       = inf.SourceSequence.Item_1.SourceIsotopeHalfLife;
        obj.RP.(actItem).source.ReferenceAirKermaRate       = inf.SourceSequence.Item_1.ReferenceAirKermaRate;
        obj.RP.(actItem).source.SourceStrengthReferenceDate = inf.SourceSequence.Item_1.SourceStrengthReferenceDate;
        obj.RP.(actItem).source.SourceType                  = inf.BrachyTreatmentType;
        
        if isfield(inf.SourceSequence.Item_1,'Private_1001_10xx_Creator')
            obj.RP.(actItem).source            = nucletronSource(inf.SourceSequence.Item_1);
            obj.RP.(actItem).source.SourceType = inf.BrachyTreatmentType;
        end
        
        
        
        
    end
end
end

function source = nucletronSource(data)

global Gy cGy deg d mm U h mCi cm seconds

source.lambda                         = data.Private_1001_1045*cGy/h/U;
source.SourceStrengthConversionFactor = data.Private_1001_1046*U/mCi;
%Raidal Dose
source.RadialDoseDistanceNumber       = data.Private_1001_1047;
source.RadialDoseDistance             = data.Private_1001_1048*mm;
source.RadialDoseValue                = data.Private_1001_1049;
%2D Anisotropy
source.AnisotropyRadialDistanceNumber  = data.Private_1001_104a;
source.AnisotropyRadialDistances       = data.Private_1001_104b*mm;
source.AnisotropyPolarAngleNumber      = data.Private_1001_104c;
source.AnisotropyPolarAngles           = data.Private_1001_104d*deg;
source.AnisotropyFunctionValue         = data.Private_1001_104e;
%1D Anisotropy
source.AnisotropyFactorRadialDistanceNumber  = data.Private_1001_104f;
source.AnisotropyFactorRadialDistance        = data.Private_1001_1050*mm;
source.AnisotropyFactorValue                 = data.Private_1001_1051;
source.AnisotropyConstant                    = data.Private_1001_1052;

%Source Strength (Air Kerma Strength)
source.SourceStrengthUnit             = data.Private_1001_1056;
source.SourceStrength                 = data.Private_1001_1055*eval(source.SourceStrengthUnit);
source.SourceStrengthImplanted        = data.ReferenceAirKermaRate*eval(source.SourceStrengthUnit);


%Source Specifivation
source.SourceIsotopeName              = data.SourceIsotopeName;
source.SourceIsotopeHalfLife          = data.SourceIsotopeHalfLife*d;
source.ActiveSourceDiameter           = data.ActiveSourceDiameter*mm;
source.ActiveSourceLength             = data.ActiveSourceLength*mm;


end



