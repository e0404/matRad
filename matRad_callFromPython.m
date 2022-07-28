function matRad_callFromPython(functionName, outputName, varargin)
%matRad_callFromPython Function that uses temporary mat file to call any function from within python

for i=1:length(varargin)
    load(varargin{i});
    var=varargin{i}(1:end-4);
    functionVars{i}=var;
end

execFunc = sprintf('%s = %s(%s);', outputName, functionName, strjoin(functionVars,','));

if functionName=='matRadGUI'
    execFunc = 'matRadGUI';
end

eval(execFunc);

if execFunc ~= 'matRadGUI'
    save(strcat(outputName,'.mat'), outputName);
end

end
