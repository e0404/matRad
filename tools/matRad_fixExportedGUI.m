fid  = fopen('matRadGUI_export.m','r');
f=fread(fid,'*char')';
fclose(fid);


f = regexprep(f,'\n''Alphamap'',\[[\s.;0-9]*?\](,\.\.\.|\))','');
f = regexprep(f,'\n''DecorationContainer'',\[[\s.;0-9]*?\](,\.\.\.|\))','');
f = regexprep(f,'\n''Layer'',''(back|middle|front)''(,\.\.\.|\))','');
f = regexprep(f,'\n''DisplayName'',blanks\(0\)(,\.\.\.|\))','');
f = regexprep(f,'\n''HelpTopicKey'',blanks\(0\)(,\.\.\.|\))','');
f = regexprep(f,'\n''DimensionNames'',{[\s\w'']*?}(,\.\.\.|\))','');
f = regexprep(f,'\n''Description'',''[\w\s]*?''(,\.\.\.|\))','');

%On/Off properties
expr = {'FontSmoothing','CLimInclude','ALimInclude','IncludeRenderer','IsContainer','IsContainer','Serializable','TransformForPrintFcnImplicitInvoke'};
expr = strjoin(expr,'|');
expr = ['\n''(' expr ')'',''(on|off)''(,\.\.\.|\))'];
f = regexprep(f,expr,'');
edit
%Mode sets
modestrings = {'ScreenPixelsPerInch','DecorationContainer','Color','FontSize','Layer','FontSmoothing','IsContainer','PickableParts','DimensionNames','Description','TransformForPrintFcnImplicitInvoke'};
modestrings = strjoin(modestrings,'|');
expr = ['\n''(' modestrings ')Mode'',''\w*?''(,\.\.\.|\))'];
f = regexprep(f,expr,'');


%Mode gets
modestrings = {'Colormap','Alphamap','Camera','DataSpace','ColorSpace','DecorationContainer','ChildContainer','XRuler','YRuler','ZRuler','AmbientLightSource','ActivePositionProperty'};
modestrings = strjoin(modestrings,'|');
expr = ['\n''(' modestrings ')Mode'',get\(0,''defaultaxes(' modestrings ')Mode''\)(,\.\.\.|\))'];
f = regexprep(f,expr,'');

%Fix whitespaces
f = regexprep(f,'(,\.\.\.)\s*?;',');');


fid  = fopen('matRadGUI_export2.m','w');
fprintf(fid,'%s',f);
fclose(fid);