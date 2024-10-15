function [data] = load_easy_file(filename)

% Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 14th October 2024
% Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html
% Open to collaboration—feel free to contact me!

% The purpose of this function is simply to load an .easy file without using pop_easy from eeglab.

% Open the file
fid = fopen(filename, 'r');

% Read the header (assuming 4 lines of header)
% header = cell(4, 1);
% for i = 1:4
%     header{i} = fgetl(fid);
% end

% Read the data
data = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'Delimiter', '\t'); % Adjust format specifiers as per your data

% Convert cell array to matrix
data = cell2mat(data);

% Close the file
fclose(fid);
end