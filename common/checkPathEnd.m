function [path] = checkPathEnd(path)

% Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 14th October 2024
% Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html
% Open to collaboration—feel free to contact me!

if path == 0; return; end
if ~endsWith(path, '/'); path = strcat(path, '/'); end

end