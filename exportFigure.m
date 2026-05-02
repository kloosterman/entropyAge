function exportFigure(fig, filename, varargin)
% EXPORTFIGURE Export figure with white background and black text.
%
% Use:
%   exportFigure(gcf, 'myfigure.png')
%   exportFigure(gcf, 'myfigure.pdf', 'Resolution', 600)
%   exportFigure(gcf, 'myfigure.png', 'FontName', 'Arial', 'FontSize', 8)
%
% Optional name-value pairs:
%   'Resolution'  : export resolution for raster formats (default = 300)
%   'FontName'    : font name to enforce (default = [])
%   'FontSize'    : font size to enforce globally (default = [])
%   'Renderer'    : renderer for figure, e.g. 'painters' (default = [])
%   'Background'  : figure/axes background color (default = 'w')

p = inputParser;
addRequired(p, 'fig', @(x) ishghandle(x, 'figure'));
addRequired(p, 'filename', @(x) ischar(x) || isstring(x));
addParameter(p, 'Resolution', 300, @isscalar);
addParameter(p, 'FontName', [], @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'FontSize', [], @(x) isempty(x) || isscalar(x));
addParameter(p, 'Renderer', [], @(x) ischar(x) || isstring(x) || isempty(x));
addParameter(p, 'Background', 'w');
parse(p, fig, filename, varargin{:});

filename   = char(p.Results.filename);
res        = p.Results.Resolution;
fontname   = p.Results.FontName;
fontsize   = p.Results.FontSize;
renderer   = p.Results.Renderer;
bg         = p.Results.Background;

if ~isempty(renderer)
    set(fig, 'Renderer', renderer);
end

% Figure background
set(fig, 'Color', bg);

% All axes, including standard axes and colorbars
ax = findall(fig, 'Type', 'axes');
for i = 1:numel(ax)
    try
        set(ax(i), 'Color', bg, ...
            'XColor', 'k', ...
            'YColor', 'k');
    end
    try
        ax(i).GridColor = [0 0 0];
        ax(i).MinorGridColor = [0 0 0];
    end
end

% Text objects
txt = findall(fig, 'Type', 'text');
set(txt, 'Color', 'k');

% Colorbars
cb = findall(fig, 'Type', 'ColorBar');
for i = 1:numel(cb)
    try
        set(cb(i), 'Color', 'k');
    end
end

% Optional font enforcement
if ~isempty(fontname)
    objs = findall(fig, '-property', 'FontName');
    set(objs, 'FontName', char(fontname));
end

if ~isempty(fontsize)
    objs = findall(fig, '-property', 'FontSize');
    set(objs, 'FontSize', fontsize);
end

drawnow;

[~,~,ext] = fileparts(lower(filename));

switch ext
    case '.png'
        exportgraphics(fig, filename, 'Resolution', res, 'BackgroundColor', bg);

    case '.jpg'
        exportgraphics(fig, filename, 'Resolution', res, 'BackgroundColor', bg);

    case '.jpeg'
        exportgraphics(fig, filename, 'Resolution', res, 'BackgroundColor', bg);

    case '.tif'
        exportgraphics(fig, filename, 'Resolution', res, 'BackgroundColor', bg);

    case '.tiff'
        exportgraphics(fig, filename, 'Resolution', res, 'BackgroundColor', bg);

    case '.pdf'
        exportgraphics(fig, filename, 'BackgroundColor', bg, 'ContentType', 'vector');

    case '.eps'
        print(fig, filename, '-depsc', sprintf('-r%d', res));

    otherwise
        error('Unsupported extension: %s', ext);
end
end