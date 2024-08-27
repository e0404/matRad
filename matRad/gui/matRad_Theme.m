classdef (Abstract) matRad_Theme
    properties (Constant, Abstract)
        name;               %Name of the Theme
        backgroundColor;    %background color in RGB [0,1]
        elementColor;       %element color in RGB [0,1]
        textColor;          %element color in RGB [0,1]
        highlightColor;     %element color in RGB [0,1]
        fontSize;           %font size
        fontWeight;         %font Weight
        fontName;           %font name
        author;             %credit for theme
        description;        %theme description
    end

    methods
        function obj = matRad_Theme()
            %As MATLAB / Octave sets the abstract properties before the constructor is called, we can validate them here
            obj.validate();
        end

        function validate(this)
            %Do not use matRad_cfg here!

            validateattributes(this.name,{'char'},{'nonempty','row'},mfilename('class'),'name');
            validateattributes(this.backgroundColor,{'double'},{'size',[1 3],'<=',1,'>=',0},mfilename('class'),'backgroundColor');
            validateattributes(this.elementColor,{'double'},{'size',[1 3],'<=',1,'>=',0},mfilename('class'),'elementColor');
            validateattributes(this.textColor,{'double'},{'size',[1 3],'<=',1,'>=',0},mfilename('class'),'textColor');
            validateattributes(this.highlightColor,{'double'},{'size',[1 3],'<=',1,'>=',0},mfilename('class'),'highlightColor');  
            validateattributes(this.fontSize,{'double'},{'scalar','positive'},mfilename('class'),'fontSize');
            validateattributes(this.fontWeight,{'char'},{'nonempty','row'},mfilename('class'),'fontWeight');
            validateattributes(this.fontName,{'char'},{'nonempty','row'},mfilename('class'),'fontName');
            validateattributes(this.author,{'char'},{'nonempty','row'},mfilename('class'),'author');
            validateattributes(this.description,{'char'},{'nonempty','row'},mfilename('class'),'description');
                       
            %Further checks on font weight
            if ~any(strcmp(this.fontWeight,{'normal','bold'}))
                error('matRad:matRad_Theme:invalidType','Invalid font weight. Must be one of ''normal'' or ''bold''');         
            end
        end

        function s = struct(this)
            s.name            = this.name;
            s.backgroundColor = this.backgroundColor;
            s.elementColor    = this.elementColor;
            s.textColor       = this.textColor;
            s.highlightColor  = this.highlightColor;
            s.fontSize        = this.fontSize;
            s.fontWeight      = this.fontWeight;
            s.fontName        = this.fontName;
            s.author          = this.author;
            s.description     = this.description;
        end
    end
end

