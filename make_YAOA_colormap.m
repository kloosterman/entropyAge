function cmap = make_YAOA_colormap(mode, n, widenMiddle, useDarkerGray)
% =========================================================
% MAKE_YAOA_COLORMAP
%
% Diverging colormap:
%   OA blue  -> gray -> YA red
%
% USAGE:
%   cmap = make_YAOA_colormap;
%   cmap = make_YAOA_colormap('lab_uniform');
%   cmap = make_YAOA_colormap('match_existing', 256, true, true);
%
% INPUTS (all optional):
%   mode           : 'lab_uniform' (default) or 'match_existing'
%   n              : number of colors (default = 256)
%   widenMiddle    : expand gray center (default = true)
%   useDarkerGray  : darker gray for contrast (default = true)
%
% REQUIREMENT (for 'match_existing'):
%   colormap_jetlightgray.mat in path
%
% =========================================================

%% Defaults
if nargin < 1 || isempty(mode), mode = 'lab_uniform'; end
if nargin < 2 || isempty(n), n = 256; end
if nargin < 3 || isempty(widenMiddle), widenMiddle = true; end
if nargin < 4 || isempty(useDarkerGray), useDarkerGray = true; end

middleGamma = 1.35;

% Anchor colors
col_Y = [1 0 0];           % YA red
col_O = [0 0.447 0.741];   % OA blue

% if useDarkerGray
%     col_mid = [0.75 0.75 0.75];
% else
%     col_mid = [0.85 0.85 0.85];
% end
% Lighter gray by default (better perceptual neutral)
if useDarkerGray
    col_mid = [0.80 0.80 0.80];   % slightly darker (optional)
else
    col_mid = [0.95 0.95 0.95];   % light gray (recommended default)
end

%% Build parameter
t = linspace(-1,1,n)';

if widenMiddle
    t = sign(t).*abs(t).^middleGamma;
end

negIdx = t <= 0;
posIdx = t > 0;

%% --------------------------------------------------------
% MODE 1: LAB uniform (recommended)
% --------------------------------------------------------
if strcmpi(mode,'lab_uniform')

    if ~(exist('rgb2lab','file') == 2 && exist('lab2rgb','file') == 2)
        error('LAB mode requires rgb2lab/lab2rgb.');
    end

    lab_O   = rgb2lab(col_O);
    lab_mid = rgb2lab(col_mid);
    lab_Y   = rgb2lab(col_Y);

    cmap_lab = zeros(n,3);

    % Blue -> gray
    tneg = (t(negIdx)+1);
    cmap_lab(negIdx,:) = (1-tneg).*lab_O + tneg.*lab_mid;

    % Gray -> red
    tpos = t(posIdx);
    cmap_lab(posIdx,:) = (1-tpos).*lab_mid + tpos.*lab_Y;

    cmap = lab2rgb(cmap_lab, 'OutputType','double');
    cmap = max(0,min(1,cmap));

%% --------------------------------------------------------
% MODE 2: match_existing (uses your jetlightgray shape)
% --------------------------------------------------------
elseif strcmpi(mode,'match_existing')

    % Load existing colormap
    S = load('colormap_jetlightgray.mat');
    fn = fieldnames(S);
    oldcmap = S.(fn{1});

    % Resample
    x_old = linspace(0,1,size(oldcmap,1));
    x_new = linspace(0,1,n);
    oldcmap_rs = interp1(x_old, oldcmap, x_new);

    % Linear RGB interpolation
    cmap = zeros(n,3);

    % Blue -> gray
    tneg = (t(negIdx)+1);
    cmap(negIdx,:) = [ ...
        col_O(1) + tneg.*(col_mid(1)-col_O(1)), ...
        col_O(2) + tneg.*(col_mid(2)-col_O(2)), ...
        col_O(3) + tneg.*(col_mid(3)-col_O(3)) ];

    % Gray -> red
    tpos = t(posIdx);
    cmap(posIdx,:) = [ ...
        col_mid(1) + tpos.*(col_Y(1)-col_mid(1)), ...
        col_mid(2) + tpos.*(col_Y(2)-col_mid(2)), ...
        col_mid(3) + tpos.*(col_Y(3)-col_mid(3)) ];

    % Transfer lightness profile (if possible)
    if exist('rgb2lab','file') == 2 && exist('lab2rgb','file') == 2
        oldLAB = rgb2lab(oldcmap_rs);
        newLAB = rgb2lab(cmap);

        Lold = oldLAB(:,1);
        mid = ceil(n/2);

        % Normalize halves
        Lleft = (Lold(1:mid)-min(Lold(1:mid))) ./ ...
                max(eps,(max(Lold(1:mid))-min(Lold(1:mid))));
        Lright = (Lold(mid:end)-min(Lold(mid:end))) ./ ...
                 max(eps,(max(Lold(mid:end))-min(Lold(mid:end))));

        Lnew = newLAB(:,1);

        L1min = min(Lnew(1:mid)); L1max = max(Lnew(1:mid));
        L2min = min(Lnew(mid:end)); L2max = max(Lnew(mid:end));

        Lnew(1:mid)   = L1min + Lleft  * (L1max-L1min);
        Lnew(mid:end) = L2min + Lright * (L2max-L2min);

        newLAB(:,1) = Lnew;
        cmap = lab2rgb(newLAB,'OutputType','double');
        cmap = max(0,min(1,cmap));
    end

else
    error('Mode must be ''lab_uniform'' or ''match_existing''.');
end
end