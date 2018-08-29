function [ pFire ] = neuralEncoder( strain,STA,NLDfun,calib )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


    strainConv = conv( [zeros(1,length(STA)-1),strain], fliplr( STA), 'valid');
    pFire = NLDfun(strainConv/calib);
    
    
    calib = max(  strainConv )
end

