function matRad_plotSOBP(ct, stf, cube1, cube2)
if exist('cube2')
    checkSize = size(cube1) == size(cube2);
    if ~all(checkSize)
        error('Dose cubes are not the same size!')
    end
end
for i = 1:length(stf)
    for k = 1:length(stf(i).ray)
        if sum(stf(i).ray(k).rayPos_bev == [0 0 0]) == 3
            middleRayIx = k;
        end
    end
    
    [alphas,l,rho,d12,ix] = matRad_siddonRayTracer(stf(i).isoCenter,ct.resolution,stf(i).sourcePoint,stf(i).ray(middleRayIx).targetPoint,{cube1});
    
    alphaMid = (alphas(1:end-1) + alphas(2:end))/2;
    physMid = alphaMid * d12;
    %physMidRel = physMid - min(physMid);
    values = cube1(ix);
    
    figure;
    plot(physMid,values)
    if exist('cube2')
        values2 = cube2(ix);
        hold on
        plot(physMid,values2)
        legend('Homogeneous','Heterogeneous')
    end
    title(['Beam ',num2str(i)])
    ylabel('RBE Dose [Gy]')
    xlabel('Distance to source')
    xlim([6500 6575])
end
end



