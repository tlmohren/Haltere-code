function [ z,E,order ] = loadComsolStrain( full_dir , which_strain )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    % import data into X (do I need if statement here? )
    if strcmp(full_dir(end-3:end),'.txt')
        X = importdata( full_dir, 't');
    else
        X = importdata( [full_dir,'.txt'] , 't');
    end

    % Check if file contains header, 
    first_line = char(X(1));
    header = all(first_line(3:7) == 'Model');
    
    % Divide string into columns 
    testline = char(X(9));
    testline_split = strsplit(testline,'  ');
    for j = 1:length(which_strain)
        matches = strfind(testline_split, which_strain{j});
        line_inds(j,:) = find(~cellfun(@isempty,matches));
    end
    
    [n_strains, n_times ] = size(line_inds);
%    
    if header
        % take out first 9 rows 
        n_points = length(X)-9;
        % initialize strain matrix 
        E =zeros(n_times,n_points,n_strains);
        
        for j = 1:n_points
            
            % split character line 
            y11 = char( X(j+9));
            y12  = strsplit( y11,'  ');

            % xyz locations
            z(1,j) = str2double(y12{1});
            z(2,j) = str2double(y12{2}); 
            z(3,j) = str2double(y12{3});

            % strain allocations
            for k = 1: n_times
                for m = 1:n_strains
                    E(k,j,m) = str2double( y12{   line_inds(m,k)    }); 
                end
            end
        end
    else
        display('Table not right format?')
    end
%     z = round(z,14);
    order = which_strain;
end