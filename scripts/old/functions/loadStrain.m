function [ xyz, E, which_strain ] = loadStrain( nameParts )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

which_strain = {'YY'};
if nameParts{4}(1:2) == 'Om'
    file = ['data' filesep strjoin(nameParts(1:3),'_') filesep strjoin(nameParts,'_') '_strain.txt'];
    X = importdata( file , 't');
    testline = char(X(9));
    testline_split = strsplit(testline,'t');
    for j = 1:length( which_strain)
        matches = strfind(testline_split, which_strain{j} );
        line_inds(j,:) = find(~cellfun(@isempty,matches));
    end
elseif nameParts{4} == 'eigenfrequency'
%     which_strain = {'YY'};
    file = ['data' filesep strjoin(nameParts(1:3),'_') filesep strjoin(nameParts,'_') '_strain.txt'];
    X = importdata( file , 't');
    testline = char(X(9));
    testline_split = strsplit(testline,'@');
    for j = 1:length( which_strain)
        matches = strfind(testline_split, which_strain{j} );
        line_inds(j,:) = find(~cellfun(@isempty,matches));
    end
end
    
header = X(1:9);
nodes = str2num(header{5}(end-7:end));
xyz = zeros(nodes,3);
    
for k = 1:nodes
    X_j  = strsplit( X{k+9} , ' ') ;
    xyz(k,:) = cellfun(@str2num, X_j(1:3) );
end
n_times = (length(X_j)-3)/6;% str2num(header{6}(end-7:end));
D = zeros(nodes,n_times, length(which_strain) );
for j = 1:length( which_strain)
    for k = 1:nodes
%         mask = ~cellfun(@isempty,matches);
%         mask = [0 , 0 , mask];
        X_j  = strsplit( X{k+9} , ' ') ;
        E(k,:,j) = cellfun(@str2num, X_j( line_inds(j,:) ) );
    end
end 