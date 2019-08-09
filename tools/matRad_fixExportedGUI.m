function matRad_fixExportedGUI(guiFile)

[filepath, name, ext] = fileparts(which(guiFile));

copyfile(guiFile,[filepath filesep name ext '.bak']);
copyfile([filepath filesep name '.fig'],[filepath filesep name '.fig.bak']);

fid  = fopen(guiFile,'r');
f=fread(fid,'*char')';
fclose(fid);


%fix strange newlines after string attributes
f = regexprep(f,'(String'',''\w+?)\r+?('')','$1$2');


f = regexprep(f,'\r?\n''Alphamap'',\[[\s.;0-9]*?\](,\.\.\.|\))','');
f = regexprep(f,'\r?\n''DecorationContainer'',\[[\s.;0-9]*?\](,\.\.\.|\))','');
f = regexprep(f,'\r?\n''Layer'',''(back|middle|front)''(,\.\.\.|\))','');
f = regexprep(f,'\r?\n''DisplayName'',blanks\(0\)(,\.\.\.|\))','');
f = regexprep(f,'\r?\n''HelpTopicKey'',blanks\(0\)(,\.\.\.|\))','');
f = regexprep(f,'\r?\n''DimensionNames'',{[\s\w'']*?}(,\.\.\.|\))','');
f = regexprep(f,'\r?\n''Description'',''[\w\s]*?''(,\.\.\.|\))','');
f = regexprep(f,'\r?\n''Tooltip'',''[\w\s]*?''(,\.\.\.|\))','');

%On/Off properties
expr = {'FontSmoothing','CLimInclude','ALimInclude','IncludeRenderer','IsContainer','IsContainer','Serializable','TransformForPrintFcnImplicitInvoke'};
expr = strjoin(expr,'|');
expr = ['\r?\n''(' expr ')'',''(on|off)''(,\.\.\.|\))'];
f = regexprep(f,expr,'');

%Mode sets
modestrings = {'ScreenPixelsPerInch','DecorationContainer','Color','FontSize','Layer','FontSmoothing','IsContainer','PickableParts','DimensionNames','Description','TransformForPrintFcnImplicitInvoke'};
modestrings = strjoin(modestrings,'|');
expr = ['\r?\n''(' modestrings ')Mode'',''\w*?''(,\.\.\.|\))'];
f = regexprep(f,expr,'');


%Mode gets
modestrings = {'Colormap','Alphamap','Camera','DataSpace','ColorSpace','FontSize','DecorationContainer','ChildContainer','XRuler','YRuler','ZRuler','AmbientLightSource','ActivePositionProperty'};
modestrings = strjoin(modestrings,'|');
expr = ['\r?\n''(' modestrings ')Mode'',get\(0,''defaultaxes(' modestrings ')Mode''\)(,\.\.\.|\))'];
f = regexprep(f,expr,'');

%Tooltip property
f = regexprep(f,'\r?\n''TooltipMode'',get\(0,''default(uipushtool|uitoggletool|uicontrol)TooltipMode''\)(,\.\.\.|\))','');

%Fix remaining whitespaces
f = regexprep(f,'(,\.\.\.)\s*?;',');');

%Fix Titles
f = regexprep(f,'''Title'',{(.*?)}','''Title'',strtrim(strjoin({$1}))');

%Octave doesn't handle ishghandle
f = regexprep(f,'ishghandle','isgraphics');

%Octave doesn't know guidemfile
f = regexprep(f,'guidemfile','%guidemfile');

%rename the MainFcn
f = regexprep(f,'gui_mainfcn',[name '_gui_mainFcn']);

%now extract the layout and mainfcn, which are added to the end of the file
%by the export
expr = '(% --- Creates and returns a handle to the GUI figure\..*?)(% --- Handles default GUIDE GUI creation and callback dispatch.*)';
out = regexp(f,expr,'tokens');
layoutFcn = out{1}{1};
guiMainFcn = out{1}{2};
%remove the functions
f = regexprep(f,expr,'');
%write the functions to files
[~,~] = mkdir([filepath filesep 'gui']);
fLayoutId = fopen([filepath  filesep 'gui' filesep name '_LayoutFcn.m'],'w');
fprintf(fLayoutId,'%s',layoutFcn);
fclose(fLayoutId);
fGuiMainFcnId = fopen([filepath filesep 'gui' filesep name '_gui_mainFcn.m'],'w');
fprintf(fGuiMainFcnId,'%s',guiMainFcn);
fclose(fGuiMainFcnId);

fid  = fopen(guiFile,'w');
fprintf(fid,'%s',f);
fclose(fid);

try 
    mat = load([filepath filesep name '.mat']);
    mat = mat.mat;
    save([filepath filesep name '.mat'],'mat','-v7');
catch
    fprintf('No .mat file was exported with GUI\n');
end

h1 = feval([name '_LayoutFcn'],'reuse');
savefig(h1,[filepath filesep name '.fig']);
close(h1);

end