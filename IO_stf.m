function status = IO_stf(path)

stf_init = load(fullfile(path,'stf_with_separate_rays.mat'), 'stf');
ray_init = load(fullfile(path,'stf_with_separate_rays.mat'), 'rays');

%%
stf = [stf_init.stf{:}];

if numel(stf) == 1
    ray_init.rays = {ray_init.rays};
end

ray = cell(size(ray_init.rays, 2), 1);

for i=1:size(ray_init.rays, 2)
    ray{i} = [ray_init.rays{i}{:}];
end

for i=1:size(ray_init.rays, 2)
    stf(i).ray = ray{i};
end

%%
save(fullfile(path,'stf.mat'), 'stf')

status = 'STF written';

end