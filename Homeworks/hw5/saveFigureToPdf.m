function saveFigureToPdf(fileName, hFig)
% saveFigureToPdf(fileName, hFig)
%
% This function saves a figure to a pdf document file
%
% INPUTS:
%   fileName = string = save the figure under this file name
%       - extension is optional (and overridden if not .pdf)
%   hFig = figure handle = optional  (default:  gcf)
% 

if nargin < 2
    hFig = gcf;
end

% Strip extension from file name
[~, fileName] = fileparts(fileName);

% Configure the figure to save in the correct format
hFig.PaperPositionMode = 'auto';
pos = hFig.PaperPosition;
hFig.PaperSize = [pos(3), pos(4)];

% Save to file
print(hFig, [fileName, '.pdf'],'-dpdf')

end

