classdef matRad_Light < matRad_Theme
    properties (Constant)
        name                = 'Light';
        backgroundColor     = [212, 211, 205]/255; %very light beige/yellow
        elementColor        = [229, 229, 229]/255; %neutral light grey
        textColor           = [ 54,  52,  46]/255; %dark yellowish-grey
        highlightColor      = [ 67,  97, 167]/255; %blue
        fontSize            = 8;
        fontWeight          = 'bold';
        fontName            = 'Helvetica';
        author              = 'Pia Stammer';
        description         = 'Default Light Theme for matRad';    
    end
end