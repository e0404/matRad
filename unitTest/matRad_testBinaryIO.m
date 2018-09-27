clear all;
load TG119.mat

%prepare metadata
%ct = evalin('base','ct');

metadata.resolution = [ct.resolution.x ct.resolution.y ct.resolution.z];
%metadata.compress = get(handles.checkbox_compress,'Value');
metadata.compress = 'false';

%Check if we have position information
% if isfield(ct,'dicomInfo')
%     if isfield(ct.dicomInfo,'ImagePositionPatient')
%        metadata.imageOrigin = ct.dicomInfo.ImagePositionPatient;       
%     end    
% end

%% Without dicom info
[metadata_post] = matRad_writeCube('./test1.nrrd',ct.cube{1},'double',metadata);

[ct_test, metadata_test] = matRad_readCube('./test1.nrrd');

assert(isequal(metadata_post.imageOrigin,metadata_test.imageOrigin),'Image Origins not equal');
assert(isequal(metadata_post.resolution,metadata_test.resolution),'Image resolutions not equal');
assert(isequal(metadata_post.axisPermutation,metadata_test.axisPermutation),'Non-matching axis permutation');
assert(isequal(metadata_post.coordinateSystem,metadata_test.coordinateSystem),'Coordinate Systems do not match');
assert(isequal(metadata_test.cubeDim,ct.cubeDim),'Wrong image dimensionality');
assert(isequal(size(ct_test),metadata_test.cubeDim),'Wrong image dimensionality');

diff = ct_test - ct.cube{1};
assert(all(diff(:) <= 1e-10), 'Significant deviations in re-imported cube!');

metadata = metadata_post;

%% With Dicom info
%Check if we have position information
if isfield(ct,'dicomInfo')
    if isfield(ct.dicomInfo,'ImagePositionPatient')
        metadata.imageOrigin = ct.dicomInfo.ImagePositionPatient;
        if ~isrow(metadata.imageOrigin)
            metadata.imageOrigin = transpose(metadata.imageOrigin);
        end
    end
end
[metadata_post] = matRad_writeCube('./test2.nrrd',ct.cube{1},'double',metadata);

[ct_test, metadata_test] = matRad_readCube('./test2.nrrd');

assert(isequal(metadata_post.imageOrigin,metadata_test.imageOrigin),'Image Origins not equal');
assert(isequal(metadata_post.resolution,metadata_test.resolution),'Image resolutions not equal');
assert(isequal(metadata_post.axisPermutation,metadata_test.axisPermutation),'Non-matching axis permutation');
assert(isequal(metadata_post.coordinateSystem,metadata_test.coordinateSystem),'Coordinate Systems do not match');
assert(isequal(metadata_test.cubeDim,ct.cubeDim),'Wrong image dimensionality');
assert(isequal(size(ct_test),metadata_test.cubeDim),'Wrong image dimensionality');

diff = ct_test - ct.cube{1};
assert(all(diff(:) <= 1e-10), 'Significant deviations in re-imported cube!');

%% With compression
metadata.compress = 'true';
[metadata_post] = matRad_writeCube('./test3.nrrd',ct.cube{1},'double',metadata);

[ct_test, metadata_test] = matRad_readCube('./test3.nrrd');

assert(isequal(metadata_post.imageOrigin,metadata_test.imageOrigin),'Image Origins not equal');
assert(isequal(metadata_post.resolution,metadata_test.resolution),'Image resolutions not equal');
assert(isequal(metadata_post.axisPermutation,metadata_test.axisPermutation),'Non-matching axis permutation');
assert(isequal(metadata_post.coordinateSystem,metadata_test.coordinateSystem),'Coordinate Systems do not match');
assert(isequal(metadata_test.cubeDim,ct.cubeDim),'Wrong image dimensionality');
assert(isequal(size(ct_test),metadata_test.cubeDim),'Wrong image dimensionality');

diff = ct_test - ct.cube{1};
assert(all(diff(:) <= 1e-10), 'Significant deviations in re-imported cube!');