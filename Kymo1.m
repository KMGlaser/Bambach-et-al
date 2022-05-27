%Matlab script for assembling kymographs
%2021-05-21
%Requires timelapse images in a folder
%First image of timelapse is displayed, in this image a rectangular region 
%for the kymograph can be defined (using the function ginputc.m). 
%Afterwards the kymograph is displayed.
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
clc
clear all
close all
%Choose folder with uigetdir:
pname=uigetdir;

%Choose which images from timelapse are used for kymograph
tmin=160;
tmax=310;
tstep=5;

%Choose min. and max. intensity displayed in the images/kymograph
C1min=0;
C1max=4000;

%-------------------------------------------------------------------------
%import timelapse images
d1=dir(pname);
cd(pname);
test_im=imread(d1(3).name);
s=size(test_im);
data=zeros(s(1),s(2),length(tmin:tstep:tmax));
count=1;
indsl=strfind(pname,'\');
tit=pname(indsl(end)+1:end);
for q=tmin:tstep:tmax
    him=imread(d1(q+3).name);
    data(:,:,count) = him;
    count=count+1;
end

%display first image
f1=figure;
im=data(:,:,1);
imshow(im,[C1min,C1max],'Border','tight','InitialMagnification',500);
axis equal
axis tight
axis off

%Choose rectangular region that is used for kymograph:
%First two points define first long edge of the rectangle, third point
%defines position of second longe edge of the rectangle
[xi1,yi1] = ginputc(1,'Color', 'g', 'LineWidth', 2, 'LineStyle', ':');
xi1=round(xi1);
yi1=round(yi1);
hold on
h1=plot(xi1,yi1,'+g','MarkerSize',10);
hold off
[xi2,yi2] = ginputc(1,'Color', 'g', 'LineWidth', 2, 'LineStyle', ':');
xi2=round(xi2);
yi2=round(yi2);
hold on
plot([xi1,xi2],[yi1,yi2],'-r')
set(h1,'Visible','Off')
hold off
[xi3,yi3] = ginputc(1,'Color', 'g', 'LineWidth', 2, 'LineStyle', ':');
v1=[xi2-xi1;yi2-yi1];
v2=[xi3-xi1;yi3-yi1];
d=abs((yi2-yi1)*xi3-(xi2-xi1)*yi3+xi2*yi1-yi2*xi1)/sqrt((yi2-yi1)^2+(xi2-xi1)^2);
if (yi1-yi2)*xi3+(xi2-xi1)*yi3-xi2*yi1+xi1*yi2>0
    alpha=pi/2;
else
    alpha=3*pi/2;
end
vp3=[cos(alpha),-sin(alpha);sin(alpha),cos(alpha)]*v1;
v3=vp3/norm(vp3)*d;
xi4=round(v3(1)+xi1);
yi4=round(v3(2)+yi1);
xi5=round(v3(1)+xi2);
yi5=round(v3(2)+yi2);
hold on
plot([xi1,xi2],[yi1,yi2],'-r')
plot([xi2,xi5],[yi2,yi5],'-r')
plot([xi5,xi4],[yi5,yi4],'-r')
plot([xi4,xi1],[yi4,yi1],'-r')
hold off

%prepare for extracting rectangular selection from timelapse images
p1=[xi1;yi1];
p2=[xi2;yi2];
p4=[xi4;yi4];
p5=[xi5;yi5];
alpha=atan((xi2-xi1)/(yi2-yi1));
s_old=size(im);
pc=s_old'/2;
rotmat=[cos(alpha),-sin(alpha);sin(alpha),cos(alpha)];
ph1=rotmat*(p1-pc);
ph2=rotmat*(p2-pc);
ph4=rotmat*(p4-pc);
ph5=rotmat*(p5-pc);
pr1=ceil(ph1+pc);
pr2=ceil(ph2+pc);
pr4=ceil(ph4+pc);
pr5=ceil(ph5+pc);
mi_x=min([pr1(1),pr2(1),pr4(1),pr5(1)]);
ma_x=max([pr1(1),pr2(1),pr4(1),pr5(1)]);
mi_y=min([pr1(2),pr2(2),pr4(2),pr5(2)]);
ma_y=max([pr1(2),pr2(2),pr4(2),pr5(2)]);
cellx=mi_y:ma_y;
if min(cellx)<=0
    disp('Warning! Selection out in rotated image!')
    cellx=cellx(cellx>0);
end
celly=mi_x:ma_x;
ss=size(data);
tmaxx=ss(3);
max_im_per_row=floor((tmax-tmin)/tstep)+1;
rows=ceil(tmaxx/max_im_per_row);
Kymo=ones(length(cellx),length(celly)*max_im_per_row,rows)*C1max;
co=0;

%extract rectangular selections from timelapse images
for q=1:tmaxx
    rotimage=imrotate(data(:,:,q),-alpha*180/pi,'bilinear','crop');
    b1=rotimage(cellx,celly);
    if co==max_im_per_row
        co=1;
    else
        co=co+1;
    end
    Kymo(1:length(cellx),((co-1)*length(celly)+1):co*length(celly),ceil(q/max_im_per_row))=b1;
end

%Plot kymograph
f2=figure;
for w=1:rows
    subplot(rows,1,w)
    imshow(Kymo(:,:,w),[C1min,C1max])
    title([tit,', tmin=',num2str(tmin),', tmax=',num2str(tmax),', tstep=',num2str(tstep)])
end


%-------------------------------------------------------------------------

function [x, y, button, ax] = ginputc(varargin)
%GINPUTC Graphical input from mouse.
%   GINPUTC behaves similarly to GINPUT, except you can customize the
%   cursor color, line width, and line style.
%
%   [X,Y] = GINPUTC(N) gets N points from the current axes and returns
%   the X- and Y-coordinates in length N vectors X and Y.  The cursor
%   can be positioned using a mouse.  Data points are entered by pressing
%   a mouse button or any key on the keyboard except carriage return,
%   which terminates the input before N points are entered.
%       Note: if there are multiple axes in the figure, use mouse clicks
%             instead of key presses. Key presses may not select the axes
%             where the cursor is.
%
%   [X,Y] = GINPUTC gathers an unlimited number of points until the return
%   key is pressed.
%
%   [X,Y] = GINPUTC(N, PARAM, VALUE) and [X,Y] = GINPUTC(PARAM, VALUE)
%   specifies additional parameters for customizing. Valid values for PARAM
%   are:
%       'Color'         : A three-element RGB vector, or one of the MATLAB
%                         predefined names, specifying the line color. See
%                         the ColorSpec reference page for more information
%                         on specifying color. Default is 'k' (black).
%       'LineWidth'     : A scalar number specifying the line width.
%                         Default is 0.5.
%       'LineStyle'     : '-', '--', '-.', ':'. Default is '-'.
%       'ShowPoints'    : TRUE or FALSE specifying whether to show the
%                         points being selected. Default is false.
%       'ConnectPoints' : TRUE or FALSE specifying whether to connect the
%                         points as they are being selected. This only
%                         applies when 'ShowPoints' is set to TRUE. Default
%                         is true.
%
%   [X,Y,BUTTON] = GINPUTC(...) returns a third result, BUTTON, that
%   contains a vector of integers specifying which mouse button was used
%   (1,2,3 from left) or ASCII numbers if a key on the keyboard was used.
%
%   [X,Y,BUTTON,AX] = GINPUTC(...) returns a fourth result, AX, that
%   contains a vector of axes handles for the data points collected.
%
%   Examples:
%       [x, y] = ginputc;
%
%       [x, y] = ginputc(5, 'Color', 'r', 'LineWidth', 3);
%
%       [x, y, button] = ginputc(1, 'LineStyle', ':');
%
%       subplot(1, 2, 1); subplot(1, 2, 2);
%       [x, y, button, ax] = ginputc;
%
%       [x, y] = ginputc('ShowPoints', true, 'ConnectPoints', true);
%
%   See also GINPUT, GTEXT, WAITFORBUTTONPRESS.

% Jiro Doke
% October 19, 2012
% Copyright 2012 The MathWorks, Inc.

try
    if verLessThan('matlab', '7.5')
        error('ginputc:Init:IncompatibleMATLAB', ...
            'GINPUTC requires MATLAB R2007b or newer');
    end
catch %#ok<CTCH>
    error('ginputc:Init:IncompatibleMATLAB', ...
        'GINPUTC requires MATLAB R2007b or newer');
end

% Check input arguments
p = inputParser();

addOptional(p, 'N', inf, @(x) validateattributes(x, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}));
addParamValue(p, 'Color', 'k', @colorValidFcn);
addParamValue(p, 'LineWidth', 0.5 , @(x) validateattributes(x, ...
    {'numeric'}, {'scalar', 'positive'}));
addParamValue(p, 'LineStyle', '-' , @(x) validatestring(x, ...
    {'-', '--', '-.', ':'}));
addParamValue(p, 'ShowPoints', false, @(x) validateattributes(x, ...
    {'logical'}, {'scalar'}));
addParamValue(p, 'ConnectPoints', true, @(x) validateattributes(x, ...
    {'logical'}, {'scalar'}));

parse(p, varargin{:});

N = p.Results.N;
color = p.Results.Color;
linewidth = p.Results.LineWidth;
linestyle = p.Results.LineStyle;
showpoints = p.Results.ShowPoints;
connectpoints = p.Results.ConnectPoints;

%--------------------------------------------------------------------------
    function tf = colorValidFcn(in)
        % This function validates the color input parameter
        
        validateattributes(in, {'char', 'double'}, {'nonempty'});
        if ischar(in)
            validatestring(in, {'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'});
        else
            assert(isequal(size(in), [1 3]) && all(in>=0 & in<=1), ...
                'ginputc:InvalidColorValues', ...
                'RGB values for "Color" must be a 1x3 vector between 0 and 1');
            % validateattributes(in, {'numeric'}, {'size', [1 3], '>=', 0, '<=', 1})
        end
        tf = true;
    end
%--------------------------------------------------------------------------

hFig = gcf;
hAx = gca;

% Save current window functions
curWBDF = get(hFig, 'WindowButtonDownFcn');
curWBMF = get(hFig, 'WindowButtonMotionFcn');
curWBUF = get(hFig, 'WindowButtonUpFcn');
curKPF  = get(hFig, 'KeyPressFcn');
curKRF  = get(hFig, 'KeyReleaseFcn');
curRF   = get(hFig, 'ResizeFcn');
try  %#ok<TRYNC> % for newer versions of MATLAB
    curWKPF = get(hFig, 'WindowKeyPressFcn');
    curWKRF = get(hFig, 'WindowKeyReleaseFcn');
end

% Save current pointer
curPointer = get(hFig, 'Pointer');
curPointerShapeCData = get(hFig, 'PointerShapeCData');

% Change window functions
set(hFig, 'WindowButtonDownFcn', @mouseClickFcn);
set(hFig, 'WindowButtonMotionFcn', @mouseMoveFcn);
set(hFig, 'WindowButtonUpFcn', '');
set(hFig, 'KeyPressFcn', @keyPressFcn);
set(hFig, 'KeyReleaseFcn', '');
set(hFig, 'ResizeFcn', @resizeFcn);
try %#ok<TRYNC> % for newer versions of MATLAB
    set(hFig, 'WindowKeyPressFcn', @keyPressFcn);
    set(hFig, 'WindowKeyReleaseFcn', '');
end

% Change actual cursor to blank
set(hFig, ...
    'Pointer', 'custom', ...
    'PointerShapeCData', nan(16, 16));

% Create an invisible axes for displaying the full crosshair cursor
hInvisibleAxes = axes(...
    'Units', 'normalized', ...
    'Position', [0 0 1 1], ...
    'XLim', [0 1], ...
    'YLim', [0 1], ...
    'HitTest', 'off', ...
    'HandleVisibility', 'off', ...
    'Visible', 'off');

% Create line object for the selected points
if showpoints
    if connectpoints
        pointsLineStyle = '-';
    else
        pointsLineStyle = 'none';
    end
    
    selectedPoints = [];
    hPoints = line(nan, nan, ...
        'Parent', hInvisibleAxes, ...
        'HandleVisibility', 'off', ...
        'HitTest', 'off', ...
        'Color', [1 0 0], ...
        'Marker', 'o', ...
        'MarkerFaceColor', [1 .7 .7], ...
        'MarkerEdgeColor', [1 0 0], ...
        'LineStyle', pointsLineStyle);
end

% Create tooltip for displaying selected points
hTooltipControl = text(0, 1, 'HIDE', ...
    'Parent', hInvisibleAxes, ...
    'HandleVisibility', 'callback', ...
    'FontName', 'FixedWidth', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', [.5 1 .5]);
hTooltip = text(0, 0, 'No points', ...
    'Parent', hInvisibleAxes, ...
    'HandleVisibility', 'off', ...
    'HitTest', 'off', ...
    'FontName', 'FixedWidth', ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', ...
    'BackgroundColor', [1 1 .5]);

% Call resizeFcn to update tooltip location
resizeFcn();

% Create full crosshair lines
hCursor = line(nan, nan, ...
    'Parent', hInvisibleAxes, ...
    'Color', color, ...
    'LineWidth', linewidth, ...
    'LineStyle', linestyle, ...
    'HandleVisibility', 'off', ...
    'HitTest', 'off');

% Prepare results
x = [];
y = [];
button = [];
ax = [];

% Wait until enter is pressed.
uiwait(hFig);


%--------------------------------------------------------------------------
    function mouseMoveFcn(varargin)
        % This function updates cursor location based on pointer location
        
        cursorPt = get(hInvisibleAxes, 'CurrentPoint');
        
        set(hCursor, ...
            'XData', [0 1 nan cursorPt(1) cursorPt(1)], ...
            'YData', [cursorPt(3) cursorPt(3) nan 0 1]);
    end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
    function mouseClickFcn(varargin)
        % This function captures mouse clicks.
        % If the tooltip control is clicked, then toggle tooltip display.
        % If anywhere else is clicked, record point.
        
        if isequal(gco, hTooltipControl)
            tooltipClickFcn();
        else
            updatePoints(get(hFig, 'SelectionType'));
        end
    end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
    function keyPressFcn(varargin)
        % This function captures key presses.
        % If "return", then exit.
        % If "delete" (or "backspace"), then delete previous point.
        % If any other key, record point.
        
        key = double(get(hFig, 'CurrentCharacter'));
        if isempty(key)
            return;
        end
        
        switch key
            case 13  % return
                exitFcn();
                
            case {8, 127}   % delete or backspace
                if ~isempty(x)
                    x(end) = [];
                    y(end) = [];
                    button(end) = [];
                    ax(end) = [];
                    
                    if showpoints
                        selectedPoints(end, :) = [];
                        set(hPoints, ...
                            'XData', selectedPoints(:, 1), ...
                            'YData', selectedPoints(:, 2));
                    end
                    
                    displayCoordinates();
                end
                
            otherwise
                updatePoints(key);
                
        end
    end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
    function updatePoints(clickType)
        % This function captures the information for the selected point
        
        hAx = gca;
        pt = get(hAx, 'CurrentPoint');
        x = [x; pt(1)];
        y = [y; pt(3)];
        ax = [ax; hAx];
        
        if ischar(clickType)   % Mouse click
            switch lower(clickType)
                case 'open'
                    clickType = 1;
                case 'normal'
                    clickType = 1;
                case 'extend'
                    clickType = 2;
                case 'alt'
                    clickType = 3;
            end
        end
        button = [button; clickType];
        
        displayCoordinates();
        
        if showpoints
            cursorPt = get(hInvisibleAxes, 'CurrentPoint');
            selectedPoints = [selectedPoints; cursorPt([1 3])];
            set(hPoints, ...
                'XData', selectedPoints(:, 1), ...
                'YData', selectedPoints(:, 2));
        end
        
        % If captured all points, exit
        if length(x) == N
            exitFcn();
        end
    end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
    function tooltipClickFcn()
        % This function toggles the display of the tooltip
        
        if strcmp(get(hTooltipControl, 'String'), 'SHOW')
            set(hTooltipControl, 'String', 'HIDE');
            set(hTooltip, 'Visible', 'on');
        else
            set(hTooltipControl, 'String', 'SHOW');
            set(hTooltip, 'Visible', 'off');
        end
    end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
    function displayCoordinates()
        % This function updates the coordinates display in the tooltip
        
        if isempty(x)
            str = 'No points';
        else
            str = sprintf('%d: %0.3f, %0.3f\n', [1:length(x); x'; y']);
            str(end) = '';
        end
        set(hTooltip, ...
            'String', str);
    end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
    function resizeFcn(varargin)
        % This function adjusts the position of tooltip when the figure is
        % resized
        
        sz = get(hTooltipControl, 'Extent');
        set(hTooltip, 'Position', [0 sz(2)]);
    end
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
    function exitFcn()
        % This function exits GINPUTC and restores previous figure settings
        
        % Restore window functions and pointer
        set(hFig, 'WindowButtonDownFcn', curWBDF);
        set(hFig, 'WindowButtonMotionFcn', curWBMF);
        set(hFig, 'WindowButtonUpFcn', curWBUF);
        set(hFig, 'KeyPressFcn', curKPF);
        set(hFig, 'KeyReleaseFcn', curKRF);
        set(hFig, 'ResizeFcn', curRF);
        set(hFig, 'Pointer', curPointer);
        set(hFig, 'PointerShapeCData', curPointerShapeCData);
        
        try %#ok<TRYNC> % for newer versions of MATLAB
            set(hFig, 'WindowKeyPressFcn', curWKPF);
            set(hFig, 'WindowKeyReleaseFcn', curWKRF);
        end
        
        % Delete invisible axes and return control
        delete(hInvisibleAxes);
        uiresume(hFig);
    end
%--------------------------------------------------------------------------

end
