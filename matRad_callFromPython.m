function matRad_callFromPython(functionName, outputName, inputPath, outputPath, varargin)
%matRad_callFromPython Function that uses temporary mat file to call any function from within python

for i=1:length(varargin)
    if contains(string(varargin{i}), string('.mat'))
        load(fullfile(inputPath, varargin{i}));
        [path, var, ext]=fileparts(varargin{i});
        functionVars{i}=var;
    else
        functionVars{i} = num2str(varargin{i});
    end
end

execFunc = sprintf('%s = %s(%s);', outputName, functionName, strjoin(functionVars,','));

%if functionName=='matRadGUI'
%    execFunc = 'matRadGUI';
%end

eval(execFunc);
save(fullfile(outputPath, [outputName '.mat']), outputName);

end
