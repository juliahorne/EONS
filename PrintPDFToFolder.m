%% Print a PDF plot and save it to a specific folder in a designated path
% Julia Horne, 2022

% This is just adding onto Colin's 'printplotpdf' function, adding in an
% extra argument that designates a folder to save the pdf file within.
% Optional argument for an out-of-path folder to save to (newpath); default
% is 'pwd' for present working directory. 

% Example usage:
%   PrintPDFToFolder(5,7,'TestFig','/TestFigures')

function PrintPDFToFolder(width,height,filename,folder)
% if the folder has not been designated, name it
if ~exist('folder','var')
    folder = 'FIGURES'; 
end

% do the basic modifications to the pdf format, given the inputted dimensions
% NOTE: units are now inches, so assume that you are printing on 8.5" x 11" 
% basic paper (recommend using width = 5.5 height = 7.5 to allow for
% caption and compensate for margins!)
set(gcf,'paperunits','inches',...
    'PaperSize',[width height],...
    'PaperPosition',[0 0 width height]);

% print the plot as a pdf with the given filename

if strcmp(folder(end),'/') && strcmp(folder(1),'/')       % don't add any extra slashes!
    FOLDER = folder;
elseif ~strcmp(folder(end),'/') && strcmp(folder(1),'/')  % make sure to add a slash after the filename
    FOLDER = [folder,'/'];
elseif strcmp(folder(end),'/') && ~strcmp(folder(1),'/')  % make sure to add a slash before the filename
    FOLDER = ['/',folder];
elseif ~strcmp(folder(end),'/') && ~strcmp(folder(1),'/') % add surrounding slashes
    FOLDER = ['/',folder,'/'];
else
    error('Incorrect folder path designation!'); 
end

% if the designated folder does not exist, then make it in this directory
stringfolder = erase(folder,'/');                         % isfolder search only works if there aren't any slashes!
if isfolder(stringfolder) == 0  
    mkdir(pwd,FOLDER); 
end

print([pwd,FOLDER,filename],'-dpdf','-r1000');             % print that PDF with a resolution of 1000 dpi - note, changing renderer to Painters will fuck with dotted/dashed lines
system(['evince ',filename,'.pdf &']);

end
