function matRad_UIInterpolation(data,dij,pln,ct,cst,optiProb)


    fInds = data.finds;
    
    
    %sort function values according to each dimension - still necessary?
    [A,I] = sort(fInds(:,1),1,'ascend');
    fIndsSorted = fInds(I,:);
    weights = data.weights(:,I);
    %idx = round(size(fInds,1)/2);
    idx = 4;
    %%
    
    slice = round(pln.propStf.isoCenter(1,3)./ct.resolution.z);
    %%

    %%
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    DosePlot = axes('Position',[.1 .5 .4 .4]);
    title('Dose Slice')
    DVHPlot = axes('Position',[.6 .5 .35 .4]);
    title('DVH');
    %
    
    %Create uipanels
    p = uipanel(f,'Position',[0.05 0.1 0.3 0.35]);
    p2 = uipanel(f,'Position',[0.35 0.3 0.15 0.15]);
    
%extract the names of the objectives and corresponding VOI
    names = {};
    for i = (1:numel(optiProb.objectives))
        names{end + 1} = optiProb.objectives{i}.name;
    end
    
    VOIs = {};
    for i = 1:size(cst,1)
        if ~isempty(cst{i,6})
            VOIs{end + 1} = cst{i,2};
        end
    end
    %Generate reference point (should lie somewhere in the middle of
    %objective 1)
    fRef = fIndsSorted(idx,:);
    wRef = weights(:,idx);
    refObj = matRad_UIData(wRef,fRef,fIndsSorted);
    
    %initial UI elements
    sliders = {};
    namesui = {};
    sliderValues = {};
    fixButtons = {};
    VOIui = {};
    
    %%Create interactive elements
    for i = 1:size(fInds,2)
        VOIui{i} =  uicontrol(p,'Style','text',...
            'Units','normalized',...
            'Position',[0.02,0.9-(i-1)*0.13, 0.15,0.1],...
            'String',VOIs{i});

        namesui{i} = uicontrol(p,'Style','text',...
            'Units','normalized',...
            'Position',[0.17,0.9-(i-1)*0.13, 0.22,0.1],...
            'String',names{i});
        
        sliderValues{i} = uicontrol(p,'Style','text',...
            'Units','normalized',...
            'Position',[0.9,0.9-(i-1)*0.13, 0.1,0.05],...
            'String',fIndsSorted(idx,i));
    
        sliders{i} = uicontrol(p,'Style','slider',...
            'Units','normalized',...
            'Min',min(fIndsSorted(:,i)),'Max',max(fIndsSorted(:,i)),...
            'Position',[0.4,0.9-(i-1)*0.13, 0.35,0.05]);
                
        
        fixButtons{i} = uicontrol(p,'Style','pushbutton',...
                'Units','normalized',...
                'Position',[0.78,0.9-(i-1)*0.13, 0.1,0.05],...
                'String','Fix');

        sliders{i}.Value = fIndsSorted(idx,i);
    
    end 
    
    DVHButton = uicontrol(p2,'Style','pushbutton',...
            'Units','normalized',...
            'Position',[0.1 0.6 0.35 0.3],...
            'String','Show DVH');
    
    ExportButton = uicontrol(p2,'Style','pushbutton',...
             'Units','normalized',...
            'Position',[0.55 0.6 0.35 0.3],...
            'String','Export to GUI');
    
    ResetConstraintButton = uicontrol(p2,'Style','pushbutton',...
             'Units','normalized',...
            'Position',[0.2 0.3 0.6 0.2],...
            'String','Reset Constraints');

    ParetoSurfaceButton = uicontrol(p2,'Style','pushbutton',...
             'Units','normalized',...
            'Position',[0.2 0.05 0.6 0.2],...
            'String','Show PS');


    %%Set callback for buttons
    for i = 1:size(fInds,2)
        set(sliders{i},'Callback',{@slider_callback,sliderValues,sliders,i,weights,dij,slice,DosePlot,refObj,optiProb,data.modcst,cst,pln,ct}) 
    end

    for i = 1:numel(fixButtons)
        set(fixButtons{i},'Callback',{@FixButton_callback,refObj,sliders,i})
    end

    set(DVHButton,'Callback',{@DVHButton_callback,DVHPlot,cst,refObj,dij,pln})
    set(ExportButton,'Callback',{@ExportButton_callback,refObj,dij})
    set(ResetConstraintButton,'Callback',{@ResetConstraintButton_callback,refObj,sliders})
    set(ParetoSurfaceButton,'Callback',{@ParetoSurfaceButton_callback,refObj})


    %initial plot of first point

    cubes = reshape(dij.physicalDose{1}*wRef,dij.doseGrid.dimensions);

    cubes = matRad_interp3(dij.doseGrid.x,dij.doseGrid.y',dij.doseGrid.z, ...
                                             cubes, ...
                                             dij.ctGrid.x,dij.ctGrid.y',dij.ctGrid.z,'linear',0);

    matRad_plotSliceWrapper(DosePlot,ct,cst,1,cubes,3,slice);
    %zoom(DosePlot,1.1);
    dvh = matRad_calcDVH(cst,cubes,'cum');
    matRad_showDVH(DVHPlot,dvh,cst,pln);

end
    % End of main file
    function DVHButton_callback(~,~,axesDVH,cst,refObj,dij,pln)
        %Shows the DVH for the current plan
        resultGUI = matRad_calcCubes(refObj.wRef,dij);
        %dvh = matRad_calcDVH(cst,doseCube,'cum');
        dvh = matRad_calcDVH(cst,resultGUI.physicalDose,'cum');

        matRad_showDVH(axesDVH,dvh,cst,pln,refObj.linestyle);
        if refObj.linestyle < 4
            refObj.linestyle = refObj.linestyle + 1;
        else
            refObj.linestyle = 1;
        end
        %
       
    end

    function ExportButton_callback(~,~,refObj,dij)
        %Export the current plan to the matRadGUi for better inspection
        resultGUI = matRad_calcCubes(refObj.wRef,dij);
        assignin('base',"resultGUI",resultGUI);
        matRadGUI;
       
    end

    
    function FixButton_callback(~,~,refObj,sliders,i)

        [lb,ub] = refObj.restrictObjective(i,sliders{i}.Value); %update refObjects bounds
        

        for i = 1:numel(sliders)
            set(sliders{i},'Min',lb(i));
            set(sliders{i},'Max',ub(i));
        end
    end


    function ResetConstraintButton_callback(~,~,refObj,sliders)
        refObj.releaseObjectiveBounds();
        
        for i = 1:numel(sliders)
            set(sliders{i},'Min',refObj.lowboundSliderInit(i));
            set(sliders{i},'Max',refObj.upboundInit(i));
        end
    end

    function ParetoSurfaceButton_callback(~,~,refObj)
    %TODO: Should only show up in low dimensions
        matRad_plotPS(refObj);
       
    end


    % Callback subfunctions to support UI actions
    function slider_callback(slider,~,textFields,sliders,idx,weights,dij,slice,DosePlot,refObj,optiProb,cstMod,cst,pln,ct)
        % Update the text for the moved slider
        %Calculate combination weights
        
        
        v = matRad_navigationProblem(refObj.fIndsAll',refObj.fRef,slider.Value,idx,refObj.upboundRed);
        %%need to check actual dimensions
        if numel(v) == 0
            %navigation algorithm didnt find a new point
            for i = 1:numel(sliders)
                set(sliders{i},'Value',refObj.fRef(i));
                set(textFields{i},'String',refObj.fRef(i));
            end

        else
        
            wnew = weights*v;%new weights
    
            fNew = matRad_objectiveFunctions(optiProb,wnew,dij,cstMod);
            fNew = optiProb.normalizeObjectives(fNew');
            refObj.fRef = fNew;
            refObj.wRef = wnew;
            
            %% should be moved to seperate function
            cubes = matRad_calcFastCubes(wnew,dij,pln);   
            %#Update the slider and text values for objective functions
        
            for i = 1:numel(sliders)
                set(sliders{i},'Value',fNew(i));
                set(textFields{i},'String',fNew(i));
            end
    
            %Plot updated plan
            matRad_plotSliceWrapper(DosePlot,ct,cst,1,cubes,3,slice);
        end
    end