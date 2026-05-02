function add_dual_colorbar_labels_figure(f, cb, left_labels, right_labels)
% Add labels on both sides of a vertical colorbar, aligned to cb.Ticks.
%
% INPUTS
%   f            : figure handle
%   cb           : colorbar handle
%   left_labels  : cell array of strings, one per tick
%   right_labels : cell array of strings, one per tick
%
% Example:
%   cb = colorbar;
%   cb.Ticks = [-3 0 3];
%   add_dual_colorbar_labels_figure(gcf, cb, {'OA','0','YA'}, {'Older','Neutral','Young'})

    drawnow

    % --- checks
    if numel(left_labels) ~= numel(cb.Ticks) || numel(right_labels) ~= numel(cb.Ticks)
        error('Number of left/right labels must match number of colorbar ticks.');
    end

    if ~strcmpi(cb.Direction, 'normal')
        warning('Colorbar Direction is reverse; labels will follow cb.Ticks positions as displayed.');
    end

    cbpos = cb.Position;   % [x y w h] in normalized figure units
    ticks = cb.Ticks(:);

    % Use limits corresponding to the colorbar scale
    lims = cb.Limits;
    if numel(lims) ~= 2 || lims(1) == lims(2)
        error('Colorbar limits are invalid.');
    end

    % Normalize tick locations to [0,1] within colorbar
    tnorm = (ticks - lims(1)) ./ (lims(2) - lims(1));

    % If reversed direction, flip positions
    if strcmpi(cb.Direction, 'reverse')
        tnorm = 1 - tnorm;
    end

    % Convert to figure-normalized y positions
    y_ticks = cbpos(2) + tnorm * cbpos(4);

    % Layout parameters
    left_w   = 0.032;
    right_w  = 0.040;
    left_gap = 0.008;
    right_gap = 0.004;
    h = 0.016;  % textbox height

    % Optional small offset if you want all labels nudged a bit
    y_offset = 0.000;

    for i = 1:numel(y_ticks)
        y0 = y_ticks(i) - h/2 + y_offset;

        annotation(f, 'textbox', ...
            [cbpos(1) - left_gap - left_w, y0, left_w, h], ...
            'String', left_labels{i}, ...
            'LineStyle', 'none', ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 6, ...
            'FitBoxToText', 'off', ...
            'Margin', 0);

        annotation(f, 'textbox', ...
            [cbpos(1) + cbpos(3) + right_gap, y0, right_w, h], ...
            'String', right_labels{i}, ...
            'LineStyle', 'none', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 6, ...
            'FitBoxToText', 'off', ...
            'Margin', 0);
    end
end