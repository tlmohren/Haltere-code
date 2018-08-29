function [ simCell] = buildSimName( baseName,parameterMat )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Haltere_Sphere_Om10_dataTestV3Clean20ms
    for j = 1:size(parameterMat,1)
       simCell{j} = ['Haltere_', parameterMat{j,1}, '_Om', num2str(parameterMat{j,2}), '_', baseName];
    end

%     simCell = baseName;
end

