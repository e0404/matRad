function helper_assignmentTest(obj,property,value)
%helper_assignmentTest Assigns property by name & value
%   This function can be used to test property set methods (e.g. for
%   exceptions
    obj.(property) = value;
end

