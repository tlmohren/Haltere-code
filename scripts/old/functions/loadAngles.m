function [ phi,theta ] = loadAngles( nameParts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%     full_dir
    
file = ['data' filesep strjoin(nameParts(1:3),'_') filesep strjoin(nameParts,'_') '_excitationAngles.txt'];
X = importdata( file , 't');
% first_line = char(X(1));
testline = char(X(9));
testline_split = strsplit(testline,'@');

which_angle = {'flapangle','theta_angle'};

header = X(1:9);
nodes = str2num(header{5}(end-7:end));
for j = 1:length(which_angle)
    matches = strfind(testline_split, which_angle{j} );
    line_inds(j,:) = find(~cellfun(@isempty,matches)) + 3;
end

X_j  = strsplit( X{10} , ' ') ;
phi = cellfun(@str2num, X_j( line_inds(1,:) ) );
theta = cellfun(@str2num, X_j( line_inds(2,:) ) );



% timeline = char(X(9));
% timeline_split = strsplit(testline,'t=');
% for j = 1:
% cellfun(@strsplit, timeline_split,{' '} );

