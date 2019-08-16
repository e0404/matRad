function matRad_fixExportedGUI(guiFile,replaceOnly)
% matRad function to postprocess figure files created from GUI. Removes
% unnecessary and compatability problem causing code and extracts the
% layout and main function storing them in separate files. Backups of the
% figure & m files are created
% 
% call
%   matRad_fixExportedGUI(guiFile,replaceOnly)
%
% input
%   guiFile:        m-file exported from guide
%   replaceOnly:    optional, default false. when set to true, only
%                   unnecessary and incompatible code is removed and no
%                   external files (i.e. layout / main function) or backups
%                   are created                 
%
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2019 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    replaceOnly = false;
end

[filepath, name, ext] = fileparts(which(guiFile));

if ~replaceOnly
    copyfile(guiFile,[filepath filesep name ext '.bak']);
    copyfile([filepath filesep name '.fig'],[filepath filesep name '.fig.bak']);
end

fid  = fopen(guiFile,'r');
f=fread(fid,'*char')';
fclose(fid);


%fix strange newlines after string attributes
f = regexprep(f,'(String'',''\w+?)\r+?('')','$1$2');
f = regexprep(f,'''String'',\[\]?','''String'',''''');


f = regexprep(f,'\r?\n''Alphamap'',\[[\s.;0-9]*?\](,\.\.\.|\))','');
f = regexprep(f,'\r?\n''DecorationContainer'',\[[\s.;0-9]*?\](,\.\.\.|\))','');
f = regexprep(f,'\r?\n''Layer'',''(back|middle|front)''(,\.\.\.|\))','');
f = regexprep(f,'\r?\n''DisplayName'',blanks\(0\)(,\.\.\.|\))','');
f = regexprep(f,'\r?\n''HelpTopicKey'',blanks\(0\)(,\.\.\.|\))','');
f = regexprep(f,'\r?\n''DimensionNames'',{[\s\w'']*?}(,\.\.\.|\))','');
f = regexprep(f,'\r?\n''Description'',''[\w\s]*?''(,\.\.\.|\))','');
f = regexprep(f,'\r?\n''Tooltip'',''[\w\s]*?''(,\.\.\.|\))','');

%uitable problems
expr = {'Units','FontUnits','RearrangeableColumns','RowStriping','ForegroundColor','Enable','HandleVisibility','FontSize','FontName','FontAngle','FontWeight'};
expr = strjoin(expr,'|');
expr = ['\r?\n''(' expr ')'',get\(0,''defaultuitable(' expr ')''\)(,\.\.\.|\))'];
f = regexprep(f,expr,'');

%On/Off properties
expr = {'FontSmoothing','CLimInclude','ALimInclude','IncludeRenderer','IsContainer','IsContainer','Serializable','TransformForPrintFcnImplicitInvoke'};
expr = strjoin(expr,'|');
expr = ['\r?\n''(' expr ')'',''(on|off)''(,\.\.\.|\))'];
f = regexprep(f,expr,'');

%Mode sets
modestrings = {'ScreenPixelsPerInch','DecorationContainer','Color','FontSize','Layer','FontSmoothing','IsContainer','PickableParts','DimensionNames','Description','TransformForPrintFcnImplicitInvoke','CurrentAxes','CurrentObject','CurrentPoint','SelectionType'};
modestrings = strjoin(modestrings,'|');
expr = ['\r?\n''(' modestrings ')Mode'',''\w*?''(,\.\.\.|\))'];
f = regexprep(f,expr,'');


%Axes Mode gets
modestrings = {'Colormap','Alphamap','Camera','DataSpace','ColorSpace','FontSize','DecorationContainer','ChildContainer','XRuler','YRuler','ZRuler','AmbientLightSource','ActivePositionProperty'};
modestrings = strjoin(modestrings,'|');
expr = ['\r?\n''(' modestrings ')Mode'',get\(0,''defaultaxes(' modestrings ')Mode''\)(,\.\.\.|\))'];
f = regexprep(f,expr,'');

%Figure Mode gets
modestrings = {'PaperSize','PaperType','PaperUnits'};
modestrings = strjoin(modestrings,'|');
expr = ['\r?\n''(' modestrings ')Mode'',get\(0,''defaultfigure(' modestrings ')Mode''\)(,\.\.\.|\))'];
f = regexprep(f,expr,'');

%Tooltip property
f = regexprep(f,'\r?\n''TooltipMode'',get\(0,''default(uipushtool|uitoggletool|uicontrol)TooltipMode''\)(,\.\.\.|\))','');

%ui element mode gets
modestrings = {'BackgroundColor','Value'};
modestrings = strjoin(modestrings,'|');
f = regexprep(f,['\r?\n''(' modestrings ')Mode'',get\(0,''default(uipushtool|uitoggletool|uicontrol)(' modestrings ')Mode''\)(,\.\.\.|\))'],'');

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

if ~replaceOnly
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
    addpath([filepath filesep 'gui']);
    fLayoutId = fopen([filepath  filesep 'gui' filesep name '_LayoutFcn.m'],'w');
    fprintf(fLayoutId,'%s',layoutFcn);
    fclose(fLayoutId);
    fGuiMainFcnId = fopen([filepath filesep 'gui' filesep name '_gui_mainFcn.m'],'w');
    fprintf(fGuiMainFcnId,'%s',guiMainFcn);
    fclose(fGuiMainFcnId);
end

fid  = fopen(guiFile,'w');
fprintf(fid,'%s',f);
fclose(fid);

if ~replaceOnly
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
end