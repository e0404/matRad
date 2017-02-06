function matRad_plotPlan3D(axesHandle,stf)
%MATRAD_PLOTPLAN3D Summary of this function goes here
%   Detailed explanation goes here

vectors = cell(numel(stf),2);

hold(axesHandle,'on');


%draws rays as lines
for fieldIx = 1:numel(stf)
    beamTarget = stf(fieldIx).isoCenter;
    beamSource = stf(fieldIx).sourcePoint + stf(fieldIx).isoCenter;
    beamVector = beamTarget - beamSource;
    
    for rayIx = 1:numel(stf(fieldIx).ray)
        ray = stf(fieldIx).ray(rayIx);
        rayTarget = ray.targetPoint + stf(fieldIx).isoCenter;
        rayVector = rayTarget - beamSource;
        
        line([beamSource(1) rayTarget(1)],[beamSource(2) rayTarget(2)],[beamSource(3) rayTarget(3)],'Parent',axesHandle,'LineStyle','-','Color',0.5*[1 1 1])
    end        
end


%{
for fieldIx = 1:numel(stf)
    rayPos_bev = cat(1,stf(fieldIx).ray.rayPos_bev);
    [~,minPosIx] = min(rayPos_bev);
    [~,maxPosIx] = max(rayPos_bev);
    
    minPosX = stf(fieldIx).ray(minPosIx(1)).rayPos + stf(fieldIx).isoCenter;
    minPosZ = stf(fieldIx).ray(minPosIx(3)).rayPos + stf(fieldIx).isoCenter;
    maxPosX = stf(fieldIx).ray(maxPosIx(1)).rayPos + stf(fieldIx).isoCenter;
    maxPosZ = stf(fieldIx).ray(maxPosIx(3)).rayPos + stf(fieldIx).isoCenter;
    
    planeX = [minPosX(1) minPosZ(1) maxPosZ(1) maxPosX(1)];
    planeY = [minPosX(2) minPosZ(2) maxPosZ(2) maxPosX(2)];
    planeZ = [minPosX(3) minPosZ(3) maxPosZ(3) maxPosX(3)];
    
    fill3(axesHandle,planeX,planeY,planeZ,'red');
end
%}
end

