function  matRadVerifyGradient(func,NumBixel)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% calculate numerical gradients

wInit = ones(NumBixel,1);

[f, g] = func(wInit);
epsilon = 1e-05;
NumRealizations = 15;

for i = 1:NumRealizations
    RandomComp = round((numel(wInit)-1).*rand(1,1));
    wDelta = wInit;
    wDelta(RandomComp) = wDelta(RandomComp) + epsilon;
    [fDelta, ~] = func(wDelta);
    numGrad = (fDelta-f)/epsilon;
    diff = ((g(RandomComp)/(numGrad))-1)*100;
    fprintf(['Component # ' num2str(RandomComp) ' - percent diff in numerical and analytical gradient = '...
        num2str(diff) '\n']);
end



end

