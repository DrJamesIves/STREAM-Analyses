function copy_easy_files(rootFolder, outputFolder)

% Author: James Ives | james.white1@bbk.ac.uk / james.ernest.ives@gmail.com
% Date: 14th October 2024
% Released under GNU GPL v3.0: https://www.gnu.org/licenses/gpl-3.0.html
% Open to collaboration—feel free to contact me!


% This function searches through all folders and subfolders of a root folder and puts them in another folder.

% Here if you don't want to use this as a function you can manually set your root and output folders
if nargin < 1
    rootFolder = 'E:\Birkbeck\STREAM\data download MW\';
    outputFolder = 'E:\Birkbeck\STREAM\Datasets\1. Raw\easy files temp\';
end

% Validate input folders
if ~isfolder(rootFolder)
    error('The specified root folder does not exist.');
end
if ~isfolder(outputFolder)
    % Create output folder if it does not exist
    mkdir(outputFolder);
end

% Recursively find and copy .easy files
% Initialize the folder stack with the root folder
folderStack = {rootFolder};

while ~isempty(folderStack)
    % Get the current folder from the stack
    currentFolder = folderStack{end};
    folderStack(end) = []; % Remove the current folder from the stack

    % Get list of .easy files in the current folder
    allFiles = dir(fullfile(currentFolder, '*.easy'));

    % Filter out files that start with a dot (.)
    easyFiles = allFiles(~startsWith({allFiles.name}, '.'));

    % Copy each .easy file to the output folder
    for i = 1:length(easyFiles)
        srcFile = fullfile(currentFolder, easyFiles(i).name);
        destFile = fullfile(outputFolder, easyFiles(i).name);
        copyfile(srcFile, destFile);
    end

    % Get list of all subfolders in the current folder
    allFolders = dir(fullfile(currentFolder, '*'));
    subfolders = allFolders([allFolders.isdir]); % Keep only directories
    subfolders = subfolders(~ismember({subfolders.name}, {'.', '..'})); % Exclude '.' and '..'

    % Add subfolders to the stack
    for i = 1:length(subfolders)
        folderStack{end+1} = fullfile(currentFolder, subfolders(i).name);
    end
end

% Display a message indicating completion
disp('All .easy files have been copied successfully.');
end