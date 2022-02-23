function [data_processed] = func_read_dpdak_cake_file(dir_data,filename,nrofscans_total)
% Reads .dat file with size n x m exported from DPDAK and returns the following n x (m+1) matrix:
% First column is 2theta, following m columns are corresponding diffraction intensities
% author:   Lukas Zauner
% contact:  lukas.zauner@tuwien.ac.at
% date:     Q4, 2021

    cur_file_dir = append(dir_data,filename);
    cur_file = table2array(readtable(cur_file_dir,'Format',repmat('%f',1,(nrofscans_total*2)),'Delimiter','\t','ReadVariableNames',false,'HeaderLines', 4));

    data_processed(:,1) = cur_file(:,1);
    data_processed(:,2:1+nrofscans_total) = cur_file(:,1+nrofscans_total:nrofscans_total*2);   
end

