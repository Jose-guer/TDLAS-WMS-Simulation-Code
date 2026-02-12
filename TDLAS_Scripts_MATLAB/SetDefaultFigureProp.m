function SetDefaultFigureProp(fontName,fontSize)
%SETDEFAULTFIGUREFONT sets the default font for figures
%   fontName is the the font name as a string and fontSize is the font
%   point size as a number

if nargin == 0
    fontName = 'Times New Roman';
    fontSize = 25;
end

% Axes
set(0,'DefaultAxesFontName', fontName);
set(0,'DefaultAxesFontSize', fontSize);
set(0,'DefaultFigureColor','w')

% Text Font
set(0,'DefaultTextFontName', fontName);
set(0,'DefaultTextFontSize', fontSize);

% LineWidth
set(0, 'DefaultLineLineWidth', 1.5);

% Interpreter
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')

% Legend
set(0, 'DefaultLegendBox', 'off')

blue = [0 0 1];
red = [1 0 0];
black = [0 0 0];
green = [76 175 80]/255;
yellow = [251 192 45]/255;
lightblue = [0 128 255]/255;
violet = [138,43,226]/255;
orange = [1 127/225 0];

newcolors = [blue; red; black; green; yellow; lightblue; violet; orange];
set(0,'DefaultAxesColorOrder',newcolors)

end

