function checkAndCreateFolders(data)

% Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 14th October 2024
% Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html
% Open to collaboration—feel free to contact me!

for p = 1:size(data, 1)
    % If data is sent through as 0 then ignore
    if data{p} == 0; continue; end
    % If path doesn't exist then make path
    if ~exist(data{p, :}, "dir"); mkdir(data{p, :}); end
end

end