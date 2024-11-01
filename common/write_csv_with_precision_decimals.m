function write_csv_with_precision_decimals(matrix_to_save, filename, decimal_places)

% ------------------------------------------------------------------------------------------------------
% Author: James Ives
% Email: james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 31st October 2024
% 
% This script was written by James Ives and is released under the GNU General Public License v3.0. 
% 
% You are free to redistribute and/or modify this script under the terms of the GNU General Public 
% License as published by the Free Software Foundation, either version 3 of the License, or (at 
% your option) any later version.
% 
% This script is provided "as-is" without any warranty; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
% details: https://www.gnu.org/licenses/gpl-3.0.html
% 
% I am happy to collaborate on any projects related to this script. 
% Feel free to contact me at the email addresses provided.
% -----------------------------------------------------------------------------------------------------

% The writematrix function in MATLAB uses the default floating-point precision when writing numeric values to a CSV file, which can lead to 
% inconsistent decimal places. MATLAB automatically trims trailing zeros, which can make numbers with fewer decimal places appear rounded.

% This script writes to csvs while enforcing a consistent number of decimal places using the fprintf function, which gives full control 
% over formatting. 

% Open the file for writing
fileID = fopen(filename, 'w');

% Format string for consistent decimal places (e.g., '%.5f' for 5 decimals)
formatSpec = repmat(['%.' num2str(decimal_places) 'f,'], 1, size(matrix_to_save, 2));
formatSpec = [formatSpec(1:end-1) '\n'];  % Replace last comma with newline

% Write each row of the matrix with consistent decimal places
for row = 1:size(matrix_to_save, 1)
    fprintf(fileID, formatSpec, matrix_to_save(row, :));
end

% Close the file
fclose(fileID);