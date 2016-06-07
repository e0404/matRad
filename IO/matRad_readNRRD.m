function [data, meta] = matRad_readNRRD(filename)
%MATRAD_READNRRD Summary of this function goes here
%   Detailed explanation goes here

hNrrdFile = fopen(filename,'r');

nrrdHeaderLine = fgetl(hNrrdFile);

%if (numel(nrrdHeaderLine) > 8)
    %error('

%fscanf(nrrdHeaderLine,'NRRD%d');




end

