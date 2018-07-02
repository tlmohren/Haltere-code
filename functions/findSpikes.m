function [ spikeInds ] = findSpikes( pFire )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%     indMax = find( pFire>=0.9);
%     dI = diff(indMax);
%     Ijump = find(dI~=1)';
% 
%     count = 1;
%     Jumps = [0,Ijump',length(indMax)];
%     for j = 1:length(Ijump)+1
%         Ival = indMax( Jumps(j)+1:Jumps(j+1));
%         spikeInds(j) = round( mean(Ival) );
%     end

    
    abs_diff=  [0, abs(diff( pFire ))];
    I_elligible =   find(  pFire >=0.8 );
    peak_start = [1,find(diff(I_elligible)>1)+1 ];
    peak_end = [find(diff(I_elligible)>1)+1 ,length(I_elligible)];
    
    for k2 = 1:length(peak_start)
        I_peak = I_elligible(peak_start(k2):peak_end(k2) );
        [V,I] = min(abs_diff(I_peak) );
        
        if ~isempty(I)
            spikeInds(k2) = I_peak(I);
        else
            spikeInds(k2) = nan;
        end
        
    end
    
end

