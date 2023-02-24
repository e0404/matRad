function status = IO_stf()

stf_init = load('C:\Users\g906c\Projects\pyRadPlan\stf_with_separate_rays.mat', 'stf');
ray_init = load('C:\Users\g906c\Projects\pyRadPlan\stf_with_separate_rays.mat', 'rays');

%%
stf = [stf_init.stf{:}];
ray = cell(size(ray_init.rays, 2), 1);

for i=1:size(ray_init.rays, 2)
    ray{i} = [ray_init.rays{i}{:}];
end

for i=1:size(ray_init.rays, 2)
    stf(i).ray = ray{i};
end

%%
save('C:\Users\g906c\Projects\pyRadPlan\stf.mat', 'stf')

status = 'STF written';

end