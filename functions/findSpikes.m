function [ spikeInds ] = findSpikes( pFire )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    indMax = find( pFire>=0.9);
    dI = diff(indMax);
    Ijump = find(dI~=1)';

    count = 1;
    Jumps = [0,Ijump',length(indMax)];
    for j = 1:length(Ijump)+1
        Ival = indMax( Jumps(j)+1:Jumps(j+1));
        spikeInds(j) = round( mean(Ival) );
    end

end

