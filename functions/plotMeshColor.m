function plotMesh( xyz, Cdata, plotStyle)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%     aa = scatter3( xyz(1,:) , xyz(2,:),  xyz(3,:));
    aa = scatter3( xyz(1,:) , xyz(2,:),  xyz(3,:),10, Cdata);
    hold on
    axis equal;
    set(aa,plotStyle{:});
        
end

