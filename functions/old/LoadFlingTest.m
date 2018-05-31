function [ eigenfreqs ] = LoadFullData( full_dir )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%     n_strains = 6;

    if strcmp(full_dir(end-3:end),'.txt')
        X = importdata( full_dir, 't');
    else

        X = importdata( [full_dir,'.txt'] , 't');
    end

    first_line = char(X(1));
    header = all(first_line(3:7) == 'Model');
    
    info = char(X(9));
    freq_line = char(X(10));
    temp = strsplit(freq_line);
    eigenfreqs = str2double(temp(4:length(temp)));
   
%     
end

