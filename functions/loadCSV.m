function [ xyz, data, header ] = loadCSV( fileName, parameters)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    % import data into X (do I need if statement here? )
% end


fid = fopen( [fileName '_data.csv'] , 'r');
for j = 1:9
    headerCell{j} = fgetl(fid);
end
fclose(fid);
header = headerCell';

display(['Loading: ' headerCell{1} ])
% headerCell'

L5  = strsplit( headerCell{5},',' );
n_points = str2num(L5{2});

L6  = strsplit( headerCell{6},',' );
L7  = strsplit( headerCell{7},'component' );

% initialization not flexible yet, figure out a way -------------
% data = zeros(n_points, n_times, length(parameters) ); 


% L9 = strsplit( headerCell{9},',' );
% matches = strfind( L9, parameters{1});
% mask = ~cellfun(@isempty,matches);
% n_times = length(nonzeros(mask) );

% for j = 
% n_parameters_in_data = 
%     mask = ~cellfun(@isempty,matches);
% n_times = ceil( str2num(L6{2}) / n_exp_repeats );

% desired_size = [n_points, n_times, n_exp_repeats ]

% data = zeros( n_points, n_times, n_exp_repeats); 
% -----------------------------------------------------------------------
L9 = strsplit( headerCell{9},',' );
for j = 1:length( parameters)
    matches = strfind( L9, parameters{j});
    mask = ~cellfun(@isempty,matches);
%     length( nonzeros( mask ))
    pre_data = csvread( [fileName '_data.csv'], 9,0);
    data( :, 1:length(nonzeros(mask)), j ) = pre_data(:, mask );
end

xyz = pre_data(:,1:3);
