function save_easy_file(filename, data)

% Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 14th October 2024
% Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html
% Open to collaboration—feel free to contact me!

% Saves the file as an easy file using the filename provided. If the filename is the same as the original it will
% overwrite the original without warning.

% Open the file
fid = fopen(filename, 'w');

% Write the data
% Assuming data is a matrix where each column corresponds to a channel
[rows, cols] = size(data);
formatSpec = repmat('%f\t', 1, cols);
formatSpec = [formatSpec(1:end-1) '\n'];
for i = 1:rows
    fprintf(fid, formatSpec, data(i, :));
end

% Close the file
fclose(fid);
end