classdef matRad_3DWidget < matRad_ViewingWidget

    % matRad_3DWidget class to generate GUI widget for 3D plan visualization 
    % 
    %
    % References
    %   -
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Copyright 2020 the matRad development team. 
    % 
    % This file is part of the matRad project. It is subject to the license 
    % terms in the LICENSE file found in the top-level directory of this 
    % distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part 
    % of the matRad project, including this file, may be copied, modified, 
    % propagated, or distributed except according to the terms contained in the 
    % LICENSE file.
    %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        viewingWidgetHandle;
    end
    
    events
        
    end
    
    methods
        %Constructor
        function this = matRad_3DWidget(viewingWidgetHandle,handleParent)
            matRad_cfg = MatRad_Config.instance();
            if nargin < 2
                handleParent = figure(...
                    'Units','normalized',...
                    'Position',[0.3 0.2 0.4 0.6],...
                    'Visible','on',...
                    'Color',matRad_cfg.gui.backgroundColor,...  'CloseRequestFcn',@(hObject,eventdata) figure1_CloseRequestFcn(this,hObject,eventdata),...
                    'IntegerHandle','off',...
                    'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
                    'MenuBar','none',...
                    'Name','MatRad 3D',...
                    'NumberTitle','off',...
                    'HandleVisibility','callback',...
                    'Tag','figure1');
                
            end
            
            this = this@matRad_ViewingWidget(handleParent);
            
            if nargin >= 1
               this.viewingWidgetHandle=viewingWidgetHandle;
              
            end
            this.lockUpdate = true;
            this.update();
               
        end
        
        function this=initialize(this)
            
        end
        
        function this=update(this,~)            
            if this.lockUpdate
                if ~isempty(this.viewingWidgetHandle) && isvalid(this.viewingWidgetHandle)
                    this.lockUpdate=false;
                    p = properties(this.viewingWidgetHandle);
                    % copy all the properties of the viewingwidget except for the widgethandle
                    for k = 1:length(p)
                        if ~strcmp(p{k},'widgetHandle') && ~strcmp(p{k},'handles') && ~strcmp(p{k},'lockUpdate')
                            try
                                this.(p{k}) = this.viewingWidgetHandle.(p{k});
                            catch
                                %warning('failed to copy property: %s', p{k});
                            end
                        end
                    end
                    this.lockUpdate=true;   
                end
                 
                this.plot3D();
            end
           
        end
        
    end
    
    methods(Access = protected)
        function this = createLayout(this)
            h88 = this.widgetHandle;
            this.createHandles();
            
        end
    end
    
    methods
           function plot3D(this)         
            
            if evalin('base','exist(''pln'')') && ...
                evalin('base','exist(''ct'')') && evalin('base','exist(''cst'')')
            
                ct  = evalin('base','ct');
                cst = evalin('base','cst');
                pln = evalin('base','pln');
                
                                
                if  evalin('base','exist(''resultGUI'')')
                    Result = evalin('base','resultGUI');
                end
                
                if evalin('base','exist(''stf'')')
                    stf = evalin('base','stf');
                else
                    stf = [];
                end
            else
                return
            end
            
            axesFig3D=axes(this.widgetHandle);
            view(axesFig3D,3);
            oldView = get(axesFig3D,'View');
            
            cla(axesFig3D);
            
            defaultFontSize = 8;
            
            %Check if we need to precompute the surface data
            if size(cst,2) < 8
                cst = matRad_computeAllVoiSurfaces(ct,cst);
                assignin('base','cst',cst);
            end
            
            set(this.widgetHandle,'Color',0.5*[1 1 1]);
            set(axesFig3D,'Color',1*[0 0 0]);
            
            %% Plot 3D structures
            hold(axesFig3D,'on');
            if this.plotContour && exist('cst','var') && exist('ct','var') %get(handles.radiobtnContour,'Value') && handles.State>0
                voiPatches = matRad_plotVois3D(axesFig3D,ct,cst,this.VOIPlotFlag,colorcube);
            end
            
            %% plot the CT slice
            if this.plotCT %get(handles.radiobtnCT,'Value')
                window = this.dispWindow{2,1}; %(2 for ct)
                ctMap = matRad_getColormap(this.ctColorMap,this.cMapSize);
                ctHandle = matRad_plotCtSlice3D(axesFig3D,ct,1,this.plane,this.slice,ctMap,window);
            end
            
            %% plot the dose slice
            if exist('Result','var')
                doseMap = matRad_getColormap(this.doseColorMap,this.cMapSize);
                doseIx  = 3;
                % if the selected display option doesn't exist then simply display
                % the first cube of the Result struct
                if ~isfield(Result,this.SelectedDisplayOption)
                    CubeNames = fieldnames(Result);
                    this.lockUpdate=false;
                    this.SelectedDisplayOption = CubeNames{1,1};
                    this.lockUpdate=true;
                end
                
                dose = Result.(this.SelectedDisplayOption);
                
                % dose colorwash
                if ~isempty(dose) && ~isvector(dose)
                    
%                     if isempty(this.dispWindow{doseIx,2})
%                         this.dispWindow{doseIx,2} = [min(dose(:)) max(dose(:))];   % set min and max dose values
%                     end
                    
                    if this.plotDose %get(handles.radiobtnDose,'Value')
                        [doseHandle,~,~] = matRad_plotDoseSlice3D(axesFig3D,ct,dose,this.plane,this.slice,this.CutOffLevel,this.doseOpacity,doseMap,this.dispWindow{doseIx,1});
                    end
                    if this.plotIsoDoseLines %get(handles.radiobtnIsoDoseLines,'Value')
                        matRad_plotIsoDoseLines3D(axesFig3D,ct,dose,this.IsoDose_Contours,this.IsoDose_Levels,this.plane,this.slice,doseMap,this.dispWindow{doseIx,1},'LineWidth',1.5);
                    end
                end
            end
            
            if this.plotPlan %get(handles.radiobtnPlan,'Value')
                matRad_plotPlan3D(axesFig3D,pln,stf);
            end
            
            %hLight = light('Parent',axesFig3D);
            %camlight(hLight,'left');
            %lighting('gouraud');
            
            xlabel(axesFig3D,'x [voxels]','FontSize',defaultFontSize)
            ylabel(axesFig3D,'y [voxels]','FontSize',defaultFontSize)
            zlabel(axesFig3D,'z [voxels]','FontSize',defaultFontSize)
            title(axesFig3D,'matRad 3D view');
            
            % set axis ratio
            ratios = [1 1 1]; %[1/ct.resolution.x 1/ct.resolution.y 1/ct.resolution.z];
            ratios = ratios([2 1 3]);
            set(axesFig3D,'DataAspectRatioMode','manual');
            set(axesFig3D,'DataAspectRatio',ratios./max(ratios));
            
            set(axesFig3D,'Ydir','reverse');
            
            set(axesFig3D,'view',oldView);
            
            %this.handles = handles;
           end
        
    end
end