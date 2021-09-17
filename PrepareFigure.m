function PrepareFigure(figureNumber, varargin)
% Adjust a figure to some predecided aspect ratio. Set the axis font,
% ticks and tick labels and the change the plot properties and save the
% figure as a single eps file or uses tikz to create an eps and tex
% file.
%
% Mandotory input arguments:
%   figureNumber: Figure number of the figure to be prepared
%
% Optional arguments:
%   The optional arguments are entered using 'PropertyName' and
%   'PropertyValue'. See below for possible properties and values.
%
% Examples of useage:
%   PrepareFigure(1);
%   PrepareFigure(gcf);
%   PrepareFigure(gcf,'PropertyName','PropertyValue')
%
%         PropertyName     PropertyValue                                                       Description
%           'fileName'    <file_name>                                   A string variable. Default: 'figure', filename of saved figure
%           'figSize'     'max'|'portrait'|[x,y,xextent,yextent]        A 4 x 1 integer vector for the figure's x, y, postion and xlength and ylength in pixels.
%                                                                       Set to 'max' to maximise figure to current screen or 'portriat' for portrait mode. Default: [90,70,1100,600].
%           'markerSize'    #1                                          An real number to set the figure marker size
%           'tikz'                                                      Call this property without a property value to save the figure as a .tikz file and in the .tex file use
%                                                                       set{\figureheight}{ycm} \set{\figurewidth}{xcm} \input(<file_name>.tikz), May need LuaLaTex if the size is too big.
%                                                                       Need matlab2tikz function (From Matlab central)
%           'pdf'                                                       Call this property without a property value to save the figure as pdf so that pdfLaTex can be used
%                                                                       Need export_fig function (From Matlab central)
%           'Beamer'                                                    Call this property without a property value to make figure for beamer with black background
%           'PowerPoint'                                                Call this property without a property value to make a figure suitable for a Power Point presentation.
%           'Word'                                                      Call this property without a property value to make a figure suitable for a Word document.
%           'Poster'                                                    Call this property without a property value to make a figure suitable for a poster
%           'noSave'                                                    Call this property to only prepare figure but not save the figure in any foramt. Used to first adjust and then add annotations
%
%                                                                       The following properties, nAxes, nPlots and nProps must be sepcified together,
%                                                                       it can be used to change properties of plots in a certain axis.
%
%           'nAxes'        [#1,#2,...]                                  A integer vector variable indiciating the axes of the plots.
%           'nPlots'       {[],[],...}                                  A cell variable, where each element is a vector of the plots
%                                                                       whose propeties will be changed, number of vectors or size of nPlots
%                                                                       must equal to number of axes or size of nAxes.
%           'nProps'       {{},{},...}                                  A cell variable, with cell elements, and each inner cell specifies the
%                                                                       standard Matlab 'property', 'value' of the plot function. Number of inner
%                                                                       cells must be equal to total number of plots specified in nPlots
%           'axisFontSize'          #1                                  To overide defalut axis value
%           'labelFontSize'         #1                                  To overide defalut axis value
%           'annotationFontSize'    #1                                  To overide defalut axis value
%           'LineWidth'             #1                                  To overide default Line width
%
% Example useage with nAxes, nPlots and nProps.
% For example if you want to change the color and marker of the 1st plot in
% axis 1 and the marker size of the 2nd plot in axes 2 (assuming
% there are two axes), you can invoke PrepareFigure as follows:
%   PrepareFigure(1,'nAxes',[1,2],'nPlots',{[1],[2]},'nProps',{{'color','red','Marker','o'},{'MarkerSize',12}}),
%
% W.D. Widanage 16-09-2012 (Lazy Sunday, listening to Metallica)


%% Set fonts for axes, labels, legends and annotations as wished

axisFont = 'Times New Roman';
axisFontWeight = 'Normal';                          % Axis font weight, 'Light', 'Normal, 'Demi','Bold'
axisFontAngle = 'Normal';                           % Axis font angle, 'Normal, 'Italic','Oblique'

idx = find(strcmp(varargin,'axisFontSize')==1);      % Axis font size in points
if isempty(idx)
    axisFontSize = 16;
else
    axisFontSize = varargin{idx+1};
end
axisFontColour = 'black';                           % Axis font colour
axisBgColour = 'white';                             % Axis background colour

labelFontWeight = 'Normal';                         % Label font weight, 'Light', 'Normal, 'Demi','Bold'
labelFontAngle = 'Normal';                          % Label font angle, 'Normal, 'Italic','Oblique'
idx = find(strcmp(varargin,'labelFontSize')==1);    % Label font size in points
if isempty(idx)
    labelFontSize = 18;
else
    labelFontSize = varargin{idx+1};
end
labelFontColour = 'black';                          % Label font colour

legendFontSize = 16;                                % Legend font size in points
legendFontColour = 'black';                         % Legend font colour

idx = find(strcmp(varargin,'annotationFontSize')==1);% Annotation font size in points
if isempty(idx)
    annotationFontSize = 18;
else
    annotationFontSize = varargin{idx+1};
end
annotationFontWeight = 'Normal';                    % Annotation font weight, 'Light', 'Normal, 'Demi','Bold'
annotationFontAngle = 'Normal';                     % Annotation font angle, 'Normal, 'Italic','Oblique'
annotationFontColour = 'black';                         % Annotation (arrow, line) colour

interpreter = 'Tex';                              % Set interpreter, Tex, Latex or None
textFontSize = 20;                                  % Set any text font size

defAxisSize = [0.13,0.21,0.775,0.650];              % Use this as default axis postion and size if only one plot, multiple plots uses matlab default
defFigSize = [90,70,1100,600];                      % Default figure postion and size
portFigSize = [90,70,800,900];                     % Portrait figure postion and size


idx = find(strcmp(varargin,'LineWidth')==1);        % Line width
if isempty(idx)
    LineWidth = 2;                                  % Use this as default line width for all plots
else
    LineWidth = varargin{idx+1};
end

idx = find(strcmp(varargin,'markerSize')==1);       % Axis font size in points
if isempty(idx)
    MarkerSize = 6;                                 % Use this as default marker size for Beamer plots
else
    MarkerSize = varargin{idx+1};
end
legendLineWidth = 1;                                % Use this as default for line width in legends
legendMarkerSize = 8;                               % Use this as default for the marker size in legends

BoxOnOff = 'off';                                   % Set to off to remove to and right edges of axes box
TickDirInOut = 'out';                               % Set ticks (mini dashes) in axes to be outwards or inwards
XGridOnOff = 'off';                                 % Set horizontal grids to on or off
YGridOnOff = 'off';                                 % Set vertical grids to on or off


% Settings for Beamer
idxBeamer = find(strcmp(varargin,'Beamer')==1);       % Check for Beamer argument
if isempty(idxBeamer)
else
    axisFont = 'Times New Roman';
    axisFontWeight = 'Normal';                          % Axis font weight, 'Light', 'Normal, 'Demi','Bold'
    axisFontAngle = 'Normal';                           % Axis font angle, 'Normal, 'Italic','Oblique'
    axisFontColour = 'white';                           % Axis font colour
    idx = find(strcmp(varargin,'axisFontSize')==1);      % Axis font size in points
    if isempty(idx)
        axisFontSize = 16;
    else
        axisFontSize = varargin{idx+1};
    end
    axisBgColour = [0.094, 0.094, 0.094];               % Axis background colour % Black-grey
    
    labelFontWeight = 'Normal';                         % Label font weight, 'Light', 'Normal, 'Demi','Bold'
    labelFontAngle = 'Normal';                          % Label font angle, 'Normal, 'Italic','Oblique'
    idx = find(strcmp(varargin,'labelFontSize')==1);    % Label font size in points
    if isempty(idx)
        labelFontSize = 16;
    else
        labelFontSize = varargin{idx+1};
    end
    labelFontColour = 'white';                          % Label font colour
    
    legendFontSize = 16;                                % Legend font size in points
    legendFontColour = 'white';                         % Legend font colour
    
    idx = find(strcmp(varargin,'annotationFontSize')==1);% Annotation font size in points
    if isempty(idx)
        annotationFontSize = 14;
    else
        annotationFontSize = varargin{idx+1};
    end
    annotationFontWeight = 'Normal';                    % Annotation font weight, 'Light', 'Normal, 'Demi','Bold'
    annotationFontAngle = 'Normal';                     % Annotation font angle, 'Normal, 'Italic','Oblique'
    annotationFontColour = 'white';                     % Annotation font and arrow box colour
    
    interpreter = 'Latex';                              % Set interpreter, Tex, Latex or None
    textFontSize = 10;                                  % Set any text font size
    
    defAxisSize = [0.13,0.21,0.775,0.650];              % Use this as default axis postion and size if only one plot, multiple plots uses matlab default
    defFigSize = [90,90,1100,650];                      % Default figure postion and size
    portFigSize = [90,70,500,900];                     % Portrait figure postion and size
    
    
    idx = find(strcmp(varargin,'LineWidth')==1);        % Line width
    if isempty(idx)
        LineWidth = 1.5;                                  % Use this as default line width for all plots
    else
        LineWidth = varargin{idx+1};
    end
    
    idx = find(strcmp(varargin,'markerSize')==1);       % Axis font size in points
    if isempty(idx)
        MarkerSize = 5;                                 % Use this as default marker size for Beamer plots
    else
        MarkerSize = varargin{idx+1};
    end
    
    legendLineWidth = 1;                                % Use this as default for line width in legends
    legendMarkerSize = 8;                               % Use this as default for the marker size in legends
end

% Settings for a Power Point slide
idxPpt = find(strcmp(varargin,'PowerPoint')==1);       % Check for PowerPoint argument
if isempty(idxPpt)
else
    axisFont = 'Bahnschrift';
    axisFontWeight = 'Normal';                      % Axis font weight, 'Light', 'Normal, 'Demi','Bold'
    axisFontAngle = 'Normal';                       % Axis font angle, 'Normal, 'Italic','Oblique'
    idx = find(strcmp(varargin,'axisFontSize')==1);      % Axis font size in points
    if isempty(idx)
        axisFontSize = 20;
    else
        axisFontSize = varargin{idx+1};
    end
    axisFontColour = 'black';                       % Axis font colour
    axisBgColour = 'white';                         % Axis background colour
    
    labelFontWeight = 'Normal';                     % Label font weight, 'Light', 'Normal, 'Demi','Bold'
    labelFontAngle = 'Normal';                      % Label font angle, 'Normal, 'Italic','Oblique
    idx = find(strcmp(varargin,'labelFontSize')==1);% Label font size in points
    if isempty(idx)
        labelFontSize = 20;
    else
        labelFontSize = varargin{idx+1};
    end
    labelFontColour = 'black';                      % Label font colour
    
    legendFontSize = 14;                            % Legend font size in points
    legendFontColour = 'black';                     % Legend font colour
    
    idx = find(strcmp(varargin,'annotationFontSize')==1);% Annotation font size in points
    if isempty(idx)
        annotationFontSize = 14;
    else
        annotationFontSize = varargin{idx+1};
    end
    annotationFontWeight = 'Normal';                % Annotation font weight, 'Light', 'Normal, 'Demi','Bold'
    annotationFontAngle = 'Normal';                 % Annotation font angle, 'Normal, 'Italic','Oblique'
    annotationFontColour = 'black';                     % Annotation (arrow, line) colour
    
    interpreter = 'tex';                            % Set interpreter, Tex, Latex or None
    
    defAxisSize = [0.13,0.15,0.775,0.770];          % Use this as default axis postion and size if only one plot, multiple plots uses matlab default
    defFigSize = [90,70,1100,600];                  % Default figure postion and size
    portFigSize = [90,70,800,900];                 % Portrait figure postion and size

    
    idx = find(strcmp(varargin,'LineWidth')==1);        % Line width
    if isempty(idx)
        LineWidth = 3;                                  % Use this as default line width for all plots
    else
        LineWidth = varargin{idx+1};
    end
    
    idx = find(strcmp(varargin,'markerSize')==1);   % Axis font size in points
    if isempty(idx)
        MarkerSize = 6;                             % Use this as default marker size for Beamer plots
    else
        MarkerSize = varargin{idx+1};
    end
    legendLineWidth = 1;                            % Use this as default for line width in legends
    legendMarkerSize = 8;                           % Use this as default for the marker size in legends
end

% Settings for a Word document
idxWord = find(strcmp(varargin,'Word')==1);         % Check for Word argument
if isempty(idxWord)
else
    axisFont = 'Times New Roman';
    axisFontWeight = 'Normal';                      % Axis font weight, 'Light', 'Normal, 'Demi','Bold'
    axisFontAngle = 'Normal';                       % Axis font angle, 'Normal, 'Italic','Oblique'
    idx = find(strcmp(varargin,'axisFontSize')==1); % Axis font size in points
    if isempty(idx)
        axisFontSize = 20;
    else
        axisFontSize = varargin{idx+1};
    end
    axisFontColour = 'black';                       % Axis font colour
    axisBgColour = 'white';                         % Axis background colour
    
    labelFontWeight = 'Normal';                     % Label font weight, 'Light', 'Normal, 'Demi','Bold'
    labelFontAngle = 'Normal';                      % Label font angle, 'Normal, 'Italic','Oblique
    idx = find(strcmp(varargin,'labelFontSize')==1);% Label font size in points
    if isempty(idx)
        labelFontSize = 25;
    else
        labelFontSize = varargin{idx+1};
    end
    labelFontColour = 'black';                      % Label font colour
    
    legendFontSize = 14;%12;%14;                    % Legend font size in points
    legendFontColour = 'black';                     % Legend font colour
    
    idx = find(strcmp(varargin,'annotationFontSize')==1);% Annotation font size in points
    if isempty(idx)
        annotationFontSize = 20;
    else
        annotationFontSize = varargin{idx+1};
    end
    annotationFontWeight = 'Normal';                % Annotation font weight, 'Light', 'Normal, 'Demi','Bold'
    annotationFontAngle = 'Normal';                 % Annotation font angle, 'Normal, 'Italic','Oblique'
    annotationFontColour = 'black';                 % Annotation (arrow, line) colour
    
    interpreter = 'Tex';                            % Set interpreter, Tex, Latex or None
    defAxisSize = [0.13,0.15,0.775,0.770];          % Use this as default axis postion and size if only one plot, multiple plots uses matlab default
    defFigSize = [90,70,1100,600];                  % Default figure postion and size
    portFigSize = [90,70,800,900];                 % Portrait figure postion and size

    
    idx = find(strcmp(varargin,'LineWidth')==1);        % Line width
    if isempty(idx)
        LineWidth = 2;                                  % Use this as default line width for all plots
    else
        LineWidth = varargin{idx+1};
    end
    idx = find(strcmp(varargin,'markerSize')==1);   % Axis font size in points
    if isempty(idx)
        MarkerSize = 6;                             % Use this as default marker size for Beamer plots
    else
        MarkerSize = varargin{idx+1};
    end
    legendLineWidth = 1;                            % Use this as default for line width in legends
    legendMarkerSize = 8;                           % Use this as default for the marker size in legends
end

% Settings for a A0 Poster document
idxPoster = find(strcmp(varargin,'Poster')==1);         % Check for Poster argument
if isempty(idxPoster)
else
    axisFont = 'Times New Roman';
    axisFontWeight = 'Normal';                      % Axis font weight, 'Light', 'Normal, 'Demi','Bold'
    axisFontAngle = 'Normal';                       % Axis font angle, 'Normal, 'Italic','Oblique'
    idx = find(strcmp(varargin,'axisFontSize')==1); % Axis font size in points
    if isempty(idx)
        axisFontSize = 20;
    else
        axisFontSize = varargin{idx+1};
    end
    axisFontColour = 'black';                       % Axis font colour
    axisBgColour = 'white';                         % Axis background colour
    
    labelFontWeight = 'Normal';                     % Label font weight, 'Light', 'Normal, 'Demi','Bold'
    labelFontAngle = 'Normal';                      % Label font angle, 'Normal, 'Italic','Oblique
    idx = find(strcmp(varargin,'labelFontSize')==1);% Label font size in points
    if isempty(idx)
        labelFontSize = 122;
    else
        labelFontSize = varargin{idx+1};
    end
    labelFontColour = 'black';                      % Label font colour
    
    legendFontSize = 18;                            % Legend font size in points
    legendFontColour = 'black';                     % Legend font colour
    
    idx = find(strcmp(varargin,'annotationFontSize')==1);% Annotation font size in points
    if isempty(idx)
        annotationFontSize = 18;
    else
        annotationFontSize = varargin{idx+1};
    end
    annotationFontWeight = 'Normal';                % Annotation font weight, 'Light', 'Normal, 'Demi','Bold'
    annotationFontAngle = 'Normal';                 % Annotation font angle, 'Normal, 'Italic','Oblique'
    annotationFontColour = 'black';                 % Annotation (arrow, line) colour
    
    interpreter = 'Tex';                            % Set interpreter, Tex, Latex or None
    defAxisSize = [0.13,0.15,0.775,0.770];          % Use this as default axis postion and size if only one plot, multiple plots uses matlab default
    defFigSize = [90,70,1100,600];                  % Default figure postion and size
    portFigSize = [90,70,800,900];                 % Portrait figure postion and size

    
    idx = find(strcmp(varargin,'LineWidth')==1);        % Line width
    if isempty(idx)
        LineWidth = 1;                                  % Use this as default line width for all plots
    else
        LineWidth = varargin{idx+1};
    end
    
    idx = find(strcmp(varargin,'markerSize')==1);   % Axis font size in points
    if isempty(idx)
        MarkerSize = 6;                             % Use this as default marker size for Beamer plots
    else
        MarkerSize = varargin{idx+1};
    end
    legendLineWidth = 1;                            % Use this as default for line width in legends
    legendMarkerSize = 8;                           % Use this as default for the marker size in legends
end
%% Get variable input arguments and set to default if uninitialised

% Check for a given file name
idx = find(strcmp(varargin,'fileName')==1);
if isempty(idx)
    fileName = 'figure';                                                % Set file name to default
else
    fileName = varargin{idx+1};
end

% Check for a given figure resolution
idx = find(strcmp(varargin,'figSize')==1);
if isempty(idx)
    set(figureNumber,'Position',defFigSize);                            % Set figure position and ascpect ration to default
else
    figSize = varargin{idx+1};
    if strcmp(figSize,'max')
        rootUnits = get(0,'Units');                                     % Get current screen units
        scrRes = get(0,'ScreenSize');                                   % Get current screen resolution
        set(figureNumber,'Position',scrRes,'Units',rootUnits)           % Set figure postion and aspect ratio to maximum screen size
    elseif strcmp(figSize,'portrait')
        set(figureNumber,'Position',portFigSize)                        % Set figure postion and aspect ratioportrait size
    else
        set(figureNumber,'Position',figSize)                            % Set figure postion and aspect ratio to give figure size
    end
end

% Check for  tikz input
idx = find(strcmp(varargin,'tikz')==1);
if isempty(idx)
    useTikz = 0;                                                     % Set useTikz to default
else
    useTikz = 1;                                       % Set useTikz value to given value
end

% Check for  pdf input
idx = find(strcmp(varargin,'pdf')==1);
if isempty(idx)
    usePDF = 0;                                                     % Set usePDF to default
else
    usePDF = 1;                                       % Set usePDF value to given value
end

% Check noSave input
idx = find(strcmp(varargin,'noSave')==1);
if isempty(idx)
    useSave = 1;                                       % Set useSave to 0
else
    useSave = 0;                                       % Set useSave to 0
end

% Graph properties
% Check for plot properties: First get the axes of the plots
idx = find(strcmp(varargin,'nAxes')==1);
if isempty(idx)
    nAxes = [];                                                         % Set nAxes to default
else
    nAxes = varargin{idx+1};                                            % Set nAxes value to given value
end

% Check for plot properties: Now get the plots
idx = find(strcmp(varargin,'nPlots')==1);
if isempty(idx)
    nPlots = [];                                                        % Set nPlots to default
else
    nPlots = varargin{idx+1};                                           % Set nPlots value to given value
end

% Check for plot properties: Now get the plot properties
idx = find(strcmp(varargin,'nProps')==1);
if isempty(idx)
    nProps = [];                                                        % Set nProps to default
else
    nProps = varargin{idx+1};                                           % Set nProps value to given value
end


%% For a given figure edit each axes, legend, annotations and all plots within an axis

% Set figure background colour
set(figureNumber,'Color',axisBgColour);


% Get number of axis in figure
hAxes = get(figureNumber,'Children');
totalAxes = numel(hAxes);

% Edit text properties. This needs to be done pror to changing X and
% Ylabels text
set(groot,'ShowHiddenHandles','on')
hText = findobj('type','text');
try                                                                 % Edit text in arrow annotations
    set(hText,...
        'FontName'   , axisFont,...
        'FontSize'   , textFontSize,...
        'FontWeight',  annotationFontWeight,...
        'FontAngle',   annotationFontAngle,...
        'Color', annotationFontColour,...
        'Interpreter',interpreter);
catch
end

set(groot,'ShowHiddenHandles','off')

if totalAxes == 1
    set(hAxes,'Position',defAxisSize)
end

% Change axes font size and axes properties
for aa = 1: totalAxes
    
    set(hAxes(aa), ...
        'FontName'  ,  axisFont    ,...
        'FontSize'  ,  axisFontSize,...
        'FontWeight',  axisFontWeight,...
        'FontAngle',   axisFontAngle,...
        'Box'       ,  BoxOnOff    , ...
        'TickDir'   ,  TickDirInOut, ...
        'TickLength', [.007 .007]  , ...
        'XMinorTick', 'off'        , ...
        'YMinorTick', 'off'        , ...
        'XGrid'     , XGridOnOff   , ...
        'YGrid'     , YGridOnOff   ,...
        'Color',  axisBgColour, ...
        'XColor', axisFontColour, ...
        'YColor', axisFontColour)
    
    
    % Set X,Y and Title
    % Change xlabel font, font size, weight and angle
    hXlabel = hAxes(aa).XLabel;
    set(hXlabel,...
        'FontName',axisFont,...
        'Interpreter',interpreter,...
        'FontSize',labelFontSize,...
        'FontWeight',labelFontWeight,...
        'FontAngle',labelFontAngle);
    
    % Change xlabel font, font size, weight and angle
    hYlabel = hAxes(aa).YLabel;
    set(hYlabel,...
        'FontName',axisFont,...
        'Interpreter',interpreter,...
        'FontSize',labelFontSize,...
        'FontWeight',labelFontWeight,...
        'FontAngle',labelFontAngle);
    
    % Change xlabel font, font size, weight and angle
    hTitle = hAxes(aa).Title;
    set(hTitle,...
        'FontName',axisFont,...
        'Interpreter',interpreter,...
        'FontSize',labelFontSize,...
        'FontWeight',labelFontWeight,...
        'FontAngle',labelFontAngle);
    
    try
        % Edit plot thickness and marker size
        hPlots = flipud(get(hAxes(aa),'Children'));
        for pp = 1:numel(hPlots)
            set(hPlots(pp), ...
                'LineWidth', LineWidth);
            set(hPlots(pp), ...
                'MarkerSize', MarkerSize);
        end
    catch
    end
end


% Edit annotation properties
set(groot,'ShowHiddenHandles','on')
hAnnotations = findobj('type','textarrow');
try                                                                 % Edit text in arrow annotations
    set(hAnnotations,...
        'FontName'   , axisFont,...
        'FontSize'   , annotationFontSize,...
        'FontWeight',  annotationFontWeight,...
        'FontAngle',   annotationFontAngle,...
        'TextColor' , annotationFontColour,...
        'Color',      annotationFontColour,...
        'Interpreter', interpreter,...
        'HorizontalAlignment','center',...
        'TextLineWidth', 1,...
        'HeadLength' , 6,...
        'HeadWidth'  , 6,...
        'HeadStyle'  ,'vback1');
catch
end



% Edit text box properties
hBox = findobj('type','textbox');
try                                                                 % Edit text in textbox annotations
    set(hBox,...
        'FontName'   , axisFont,...
        'FontSize'   , annotationFontSize,...
        'FontWeight',  annotationFontWeight,...
        'FontAngle',   annotationFontAngle,...
        'Color' , annotationFontColour,...
        'EdgeColor', annotationFontColour,...
        'Interpreter', interpreter,...
        'HorizontalAlignment','center',...
        'LineStyle','none');
catch
end

% Edit text box properties
hText = findobj('type','text');
% try                                                                 % Edit text in textbox annotations
%     set(hText,...
%         'FontName'   , axisFont,...
%         'FontSize'   , annotationFontSize,...
%         'FontWeight',  annotationFontWeight,...
%         'FontAngle',   annotationFontAngle,...
%         'Color' , annotationFontColour,...
%         'EdgeColor', annotationFontColour,...
%         'Interpreter', interpreter,...
%         'HorizontalAlignment','center',...
%         'LineStyle','none');
% catch
% end

% Edit double end arrow properties
hBox = findobj('type','DoubleEndArrow');
set(hBox,...
    'Color' , annotationFontColour);

set(groot,'ShowHiddenHandles','off')

% Get legend handle
legendHandles = findobj(figureNumber,'Tag','legend');                   % Get the handles of the legend axes

% Edit legend properties and remove legend handles from axesHandles to get the xy axes handles
for ii = 1: length(legendHandles)                                       % Loop over each legend
    idx = axesHandles == legendHandles(ii);
    axesHandles(idx) = [];
    
    % Change legend properties
    set(legendHandles(ii),...
        'FontName', axisFont,...                                        % Change legend font
        'FontSize', legendFontSize,...                                  % Change legend font size
        'Interpreter', interpreter);                                    % Set legend interpreter to Latex
    
    
    % Get the string handles of within each legend
    legendStringHandle_temp = get(legendHandles(ii),'Children');
    lengthTemp = length(legendStringHandle_temp);
    
    for jj = 1: lengthTemp/3                                             % Every thrid handle is the string handle, the other two are line handles
        legendStringHandle = legendStringHandle_temp(jj*3);
        idxText = textHandles == legendStringHandle;                     % Find the corresponding text handle in the string handle and delete
        textHandles(idxText) = [];
        set(legendStringHandle_temp(jj*3-1),'LineWidth',legendLineWidth);   % Set legend line width
        set(legendStringHandle_temp(jj*3-2),'MarkerSize',legendMarkerSize); % Set marker size
    end
    
end


% Save file as an eps
if useSave
    if useTikz
        matlab2tikz('figurehandle',figureNumber,[fileName,'.tex'], 'height', '\figureheight', 'width', '\figurewidth','showInfo',false);
    elseif usePDF
        export_fig(figureNumber, fileName, '-pdf', '-transparent') % Save as a pdf
    elseif ~isempty(idxBeamer)
        export_fig(figureNumber, fileName, '-pdf', '-transparent') % Save as a pdf
    elseif isempty(idxPpt) && isempty(idxWord)      % An eps figure is not created when prepared for Power Point or Word
        set(figureNumber,'PaperPositionMode','auto')
        print(figureNumber,'-depsc2',fileName)      % Save as an eps
        disp(['Figure is saved with filename: ',fileName,'.eps']);
    end
    savefig(figureNumber,fileName);                               % Also save the .fig file
end

