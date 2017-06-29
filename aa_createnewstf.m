vec = [];
for i = 1:stf.numOfRays
    vec = [vec stf.ray(i).energy];
end

univec = unique(vec);

for i = 1:length(univec)
    strct.rayPerEnergy{i,1} = univec(i);
    strct.rayPerEnergy{i,2} = [];
    for j = 1:stf.numOfRays
        for k =  1:length(stf.ray(j).energy)
            if any(stf.ray(j).energy(k)==univec(i)) 
                strct.rayPerEnergy{i,2} = [strct.rayPerEnergy{i,2}; stf.ray(j).rayPos_bev];
            end
        end
    end
end



