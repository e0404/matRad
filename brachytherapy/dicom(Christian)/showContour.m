function showContour(obj,rsItem,structName)
%SHOWCONTOUR Summary of this function goes here
%   Detailed explanation goes here

clf
hold on

if nargin<3
structName = fieldnames(obj.RS.(rsItem).structures);
else
    structName = cellstr(structName);
end

for is=1:length(structName)
for im=1:length(obj.RS.(rsItem).structures.(structName{is}).outerContourMap)
    noContours = numel(obj.RS.(rsItem).structures.(structName{is}).outerContourMap(im).contourMap);
    for ic=1:noContours
        contours = obj.RS.(rsItem).structures.(structName{is}).outerContourMap(im).contourMap{ic};
        plot3(contours(:,1),contours(:,2),contours(:,3),'Color',obj.RS.(rsItem).structures.(structName{is}).color)
    end
end
end
legend(structName,'Interpreter','none')
view(45,45)
grid on
xlabel('x')
ylabel('y')
zlabel('y')



end

