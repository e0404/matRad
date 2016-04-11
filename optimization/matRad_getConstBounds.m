function [cl,cu] = matRad_getConstBounds(constraint,param)


if isequal(constraint.type, 'max dose constraint') 

    cl = [cl;-inf];
    cu = [cu;param];

elseif isequal(constraint.type, 'min dose constraint') 

    cl = [cl;param];
    cu = [cu;inf];

elseif isequal(constraint.type, 'min mean dose constraint') 

    cl = [cl;param];
    cu = [cu;inf];

elseif isequal(constraint.type, 'max mean dose constraint') 

    cl = [cl;-inf];
    cu = [cu;param];

elseif isequal(constraint.type, 'min max mean dose constraint') 

    cl = [cl;param(1)];
    cu = [cu;param(2)];

elseif isequal(constraint.type, 'min EUD constraint') 

    cl = [cl;param];
    cu = [cu;inf];

elseif isequal(constraint.type, 'max EUD constraint') 

    cl = [cl;-inf];
    cu = [cu;param];

elseif isequal(constraint.type, 'min max EUD constraint') 

    cl = [cl;param(1)];
    cu = [cu;param(2)];

elseif isequal(constraint.type, 'exact DVH constraint')

    cl = [cl;constraint.volume/100];
    cu = [cu;constraint.volume/100];

elseif isequal(constraint.type, 'max DVH constraint') 

    cl = [cl;-inf];
    cu = [cu;constraint.volume/100];

    % alternative constraint calculation 1/4 %                
    % cl = [cl;-inf];
    % cu = [cu;0];
    % alternative constraint calculation 1/4 %

elseif isequal(constraint.type, 'min DVH constraint') 

    cl = [cl;constraint.volume/100];
    cu = [cu;inf];

    % alternative constraint calculation 2/4 %                
    % cl = [cl;-inf];
    % cu = [cu;0];
    % alternative constraint calculation 2/4 %

end % constraint switch
