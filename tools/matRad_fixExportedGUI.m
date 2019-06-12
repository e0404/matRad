fid  = fopen('matRadGUI_export.m','r');
f=fread(fid,'*char')';
fclose(fid);

if ispc
    nl = '\r\n';%newline;
else
    nl = '\n';
end

f = regexprep(f,[nl '''Alphamap'',\[[\s.;0-9]*?\](,\.\.\.|\))'],'');
f = regexprep(f,[nl '''DecorationContainer'',\[[\s.;0-9]*?\](,\.\.\.|\))'],'');
f = regexprep(f,[nl '''Layer'',''(back|middle|front)''(,\.\.\.|\))'],'');
f = regexprep(f,[nl '''DisplayName'',blanks\(0\)(,\.\.\.|\))'],'');
f = regexprep(f,[nl '''HelpTopicKey'',blanks\(0\)(,\.\.\.|\))'],'');
f = regexprep(f,[nl '''DimensionNames'',{[\s\w'']*?}(,\.\.\.|\))'],'');
f = regexprep(f,[nl '''Description'',''[\w\s]*?''(,\.\.\.|\))'],'');

%On/Off properties
expr = {'FontSmoothing','CLimInclude','ALimInclude','IncludeRenderer','IsContainer','IsContainer','Serializable','TransformForPrintFcnImplicitInvoke'};
expr = strjoin(expr,'|');
expr = [nl '''(' expr ')'',''(on|off)''(,\.\.\.|\))'];
f = regexprep(f,expr,'');

%Mode sets
modestrings = {'ScreenPixelsPerInch','DecorationContainer','Color','FontSize','Layer','FontSmoothing','IsContainer','PickableParts','DimensionNames','Description','TransformForPrintFcnImplicitInvoke'};
modestrings = strjoin(modestrings,'|');
expr = [nl '''(' modestrings ')Mode'',''\w*?''(,\.\.\.|\))'];
f = regexprep(f,expr,'');


%Mode gets
modestrings = {'Colormap','Alphamap','Camera','DataSpace','ColorSpace','DecorationContainer','ChildContainer','XRuler','YRuler','ZRuler','AmbientLightSource','ActivePositionProperty'};
modestrings = strjoin(modestrings,'|');
expr = [nl '''(' modestrings ')Mode'',get\(0,''defaultaxes(' modestrings ')Mode''\)(,\.\.\.|\))'];
f = regexprep(f,expr,'');

%Fix whitespaces
f = regexprep(f,'(,\.\.\.)\s*?;',');');

%Disable compability error
f = regexprep(f,'(matRad GUI not available for)(.*?)(return;)','matRad GUI not fully supported for$2');

%Fix Javax
f = regexprep(f,'(javax\.swing\.UIManager\.setLookAndFeel\(lf\);)','try, $1, catch, fprintf(''javax not supported\\n''), end');

%Fix Movegui
f = regexprep(f,'(movegui\(gui_hFigure,''onscreen''\);)','try, $1, catch, fprintf(''movegui not supported\\n''), end');

%Fix OpenGL
f = regexprep(f,'(opengl software)','try, $1, catch, fprintf(''opengl not supported\\n''), end');

%Fix Logos
f = regexprep(f,'(matrad\_logo\.png|DKFZ\_Logo\.png)','gui/$1');

%Fix cell array titles
f = regexprep(f,'(''Title'',)({.*?})','$1strtrim(strjoin($2))');

fid  = fopen('matRadGUI_export2.m','w');
fprintf(fid,'%s',f);
fclose(fid);