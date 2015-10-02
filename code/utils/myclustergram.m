function varargout = myclustergram(data,varargin)
%CLUSTERGRAM creates a dendrogram and heat map on the same figure.
%
%   CLUSTERGRAM(DATA) creates a dendrogram and heat map from DATA using
%   hierarchical clustering with Euclidean distance metric and average
%   linkage used to generate the hierarchical tree. The clustering is
%   performed on the rows of DATA. The rows of DATA are typically genes and
%   the column are the results from different microarrays.  To cluster the
%   columns instead of the rows, transpose the data using the ' operator.
%
%   CLUSTERGRAM(...,'ROWLABELS',ROWLABELS) uses the contents of cell array
%   of strings, or a numeric array, ROWLABELS as labels for the rows in
%   DATA.
%
%   CLUSTERGRAM(...,'COLUMNLABELS',COLUMNLABELS) uses the contents of cell
%   array of strings, or a numeric array, COLUMNLABELS as labels for the
%   columns in DATA.
%
%   CLUSTERGRAM(...,'PDIST',DISTANCE) allows you to set the distance metric
%   used by PDIST, the function used to calculate pairwise distance between
%   observations. If the distance metric requires extra arguments then
%   these should be passed as a cell array, for example to use the
%   Minkowski distance with exponent P you would use {'minkowski', P}. See
%   the help for PDIST for more details of the available options. The
%   default distance metric for CLUSTERGRAM is 'Euclidean'.
%
%   CLUSTERGRAM(...,'LINKAGE',METHOD) allows you to select the linkage
%   method used by LINKAGE, the function used to create the hierarchical
%   cluster tree. See the help for LINKAGE for more details of the
%   available options. The default linkage method used by CLUSTERGRAM is
%   'average' linkage.
%
%   CLUSTERGRAM(...,'DENDROGRAM',DENDARGS) allows you to pass arguments to
%   the DENDROGRAM, the function used to create the dendrogram. DENDARGS
%   should be a cell arrays of parameter/value pairs that can be passed to
%   DENDROGRAM. See the help for DENDROGRAM for more details of the
%   available options.
%
%   CLUSTERGRAM(...,'OPTIMALLEAFORDER',FALSE) disables the optimal
%   leaf-ordering calculation. When working with large data sets,
%   calculating the optimal leaf ordering can be very time consuming and
%   uses a large amount of memory. This option is disabled by default when
%   the number of rows or columns is greater than 1000. Set the value to
%   TRUE to override this default.
%
%   CLUSTERGRAM(...,'COLORMAP',CMAP) allows you to specify the colormap
%   that will be used for the figure containing the clustergram. This will
%   control the colors used to display the heat map. It can be the name or
%   function  handle of a function that returns a colormap, or an M by 3
%   array  containing RGB values. The default is REDGREENCMAP.
%
%   CLUSTERGRAM(...,'SYMMETRICRANGE',FALSE) disables the default behavior
%   of forcing the color scale of the heat map to be symmetric about zero.
%
%   CLUSTERGRAM(...,'DIMENSION',DIM) allows you to specify whether to
%   create a  one-dimensional or two-dimensional clustergram. The options
%   are 1 or 2. The default value is 1. The one-dimensional clustergram
%   will cluster the  rows of the data. The two-dimensional clustergram
%   creates the one-dimensional  clustergram, and then clusters the columns
%   of the row-clustered data.
%
%   CLUSTERGRAM(...,'RATIO',R) allows you to specify the ratio of the space
%   that the dendrogram(s) will use, relative to the size of the heat map,
%   in the X and Y directions. If R is a single scalar value, it is used as
%   the ratio for both directions. If R is a two-element vector, the first
%   element is used for the X ratio, and the second element is used for the
%   Y ratio. The Y ratio is ignored for one-dimensional clustergrams. The
%   default ratio is 1/5.
%
%   Hold the mouse button down over the image to see the exact values at a
%   particular point.
%
%   Examples:
%
%       load filteredyeastdata;
%       clustergram(yeastvalues);
%
%       % Add some labels.
%       clustergram(yeastvalues,'ROWLABELS',genes,'COLUMNLABELS',times);
%
%       % Change the clustering parameters.
%       clustergram(yeastvalues,'PDIST','correlation','LINKAGE','complete');
%
%       % Change the dendrogram color parameter.
%       clustergram(yeastvalues,'ROWLABELS',genes,'DENDROGRAM',{'color',5});
%
%   See also DENDROGRAM, LINKAGE, MAPCAPLOT, PDIST, REDGREENCMAP.

%   Reference:
%   Bar-Joseph Z, Gifford DK, Jaakkola TS. Fast optimal leaf ordering for
%   hierarchical clustering. Bioinformatics. 2001;17 Suppl 1:S22-9. PMID:
%   11472989

% Copyright 2003-2006 The MathWorks, Inc.
% $Revision: 1.14.6.18 $   $Date: 2006/12/12 23:56:12 $

if  nargin < 1
    % get some data -- makeRandomData creates somewhat structured data that
    % looks pretty when clustered
    [data, rowLabels, colLabels ] = makeRandomData;
else
    colLabels = cell(size(data,2),1);
    rowLabels = cell(size(data,1),1);
end

% defaults
clustStruct.data = data;
clustStruct.rowLabels = rowLabels;
clustStruct.colLabels = colLabels;
clustStruct.pdistargs = {'euc'};
clustStruct.linkageargs = {'average'};
%clustStruct.dendroargs = {};
clustStruct.dendroargs = {'colorthreshold',0.7};
clustStruct.colormap = redgreencmap(256);
clustStruct.ratio = [1/5 1/5];
clustStruct.onedimension = true;
clustStruct.fullwidth = .72;
clustStruct.fullheight = .85;
clustStruct.zeroCenter = true;
clustStruct.reuse = false;
clustStruct.optLeafOrder = max(size(data)) <= 1000;
clustStruct.saturation = inf;

if nargin > 1
    if rem(nargin,2) == 0
        error('Bioinfo:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'rowlabels','columnlabels','pdist','linkage','dendrogram',...
        'colormap','dimension','ratio','symmetricrange','reuse',...
        'optimalleaforder','saturation'};
    for j=1:2:nargin-1
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error('Bioinfo:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:AmbigousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1  % rowlabels
                    clustStruct.rowLabels = pval;
                    if isnumeric(rowLabels)
                        clustStruct.rowLabels = cellstr(num2str(clustStruct.rowLabels(:)));
                    end
                case 2  % column labels
                    clustStruct.colLabels = pval;
                    if isnumeric(clustStruct.colLabels)
                        clustStruct.colLabels = cellstr(num2str(clustStruct.colLabels(:)));
                    end
                case 3  % pdist args
                    if iscell(pval)
                        clustStruct.pdistargs = pval;
                    else
                        clustStruct.pdistargs = {pval};
                    end
                case 4  % linkage args
                    if iscell(pval)
                        clustStruct.linkageargs = pval;
                    else
                        clustStruct.linkageargs = {pval};
                    end
                case 5  % dendrogram args
                    if iscell(pval)
                        clustStruct.dendroargs = pval;
                    else
                        clustStruct.dendroargs = {pval};
                    end
                case 6 % colormap

                    validcmap = true;
                    if ischar(pval)
                        cmap = eval(pval);
                    elseif isa(pval,'function_handle')
                        cmap = feval(pval);
                    elseif isnumeric(pval)
                        cmap = pval;
                    else
                        validcmap = false;
                    end

                    if ~validcmap || ndims(cmap) ~= 2 || size(cmap,2) ~= 3
                        error('Bioinfo:InvalidClustergramColormap','Invalid colormap')
                    end

                    clustStruct.colormap = cmap;
                case 7 % dimension args
                    if pval == 2
                        clustStruct.onedimension = false;
                    elseif pval == 1
                        clustStruct.onedimension = true;
                    else
                        error('Bioinfo:InvalidClustergramDimension','Invalid dimension specified: %d',pval);
                    end
                case 8 % ratio args
                    numpval = numel(pval);
                    if ~isa(pval,'double')
                        error('Bioinfo:InvalidClustergramRatioClass','Ratio must be a double.');
                    elseif numpval ~= 1 && numpval ~= 2
                        error('Bioinfo:InvalidClustergramRatioSize','Ratio must have 1 or 2 elements.');
                    elseif ~all(pval > 0 & pval < 1);
                        error('Bioinfo:InvalidClustergramRatioValue','Ratio must be between 0 and 1.')
                    end
                    if numpval == 1
                        clustStruct.ratio = [pval pval];
                    else
                        clustStruct.ratio = pval;
                    end
                case 9 % symmetric range -- whether or not to center color range about 0
                    zeroCenterFlag = opttf(pval);
                    if ~isa(pval,'logical') || numel(pval) ~= 1
                        error('Bioinfo:InputOptionNotLogicalScalar','%s must be a logical scalar, true or false.',...
                            upper(char(okargs(k))));
                    end
                    clustStruct.zeroCenter = zeroCenterFlag;
                case 10 % reuse - whether to reuse the same figure
                    if ~isa(pval,'logical') || numel(pval) ~= 1
                        error('Bioinfo:InputOptionNotLogicalScalar','%s must be a logical scalar, true or false.',...
                            upper(char(okargs(k))));
                    end
                    clustStruct.reuse = pval;
                case 11 % optimal leaf ordering
                    if ~isa(pval,'logical') || numel(pval) ~= 1
                        error('Bioinfo:InputOptionNotLogicalScalar','%s must be a logical scalar, true or false.',...
                            upper(char(okargs(k))));
                    end
                    clustStruct.optLeafOrder = pval;
                case 12 % saturation value
                    if ~isnumeric(pval) || numel(pval) ~= 1 || pval <= 0
                        error('Bioinfo:clustergram:saturationNotNumeric',...
                            'Saturation must be a positive numeric scalar value.');
                    end
                    clustStruct.saturation = pval;
            end
        end
    end
end


hFig = [];
hFigs = findall(0,'Type','figure','Tag','ClustergramFigure');
if clustStruct.reuse && ~isempty(hFigs)
    % use GCF to find a figure handle to use

    if isempty(hFigs)
        % no figures found, use GCF, which will
        % create a figure if necessary
        hFig = gcf;
        set(hFig,'Tag','ClustergramFigure');
    else
        % loop over the figures that were found, find one to reuse
        for n = 1:length(hFigs)
            tmpclustStruct = getappdata(hFigs(n),'ClustergramStructure');
            if tmpclustStruct.reuse
                hFig = hFigs(n);
                break
            end
        end
        if isempty(hFig)
            hFig = gcf;
            set(gcf,'Tag','ClustergramFigure');
        end
    end

    % delete the axes in the figure
    delete(findall(hFig,'type','axes'));
else
    hFig = gcf;
    delete(findall(hFig,'type','axes'))
    set(hFig,'Tag','ClustergramFigure');
end

% bring figure forward
figure(hFig);

if clustStruct.onedimension
    [hAxes,hImage,hDen] = create1DClustergram(clustStruct,hFig);
else
    [hAxes,hImage,hDen] = create2DClustergram(clustStruct,hFig);
end

% set the HandleVisibility for the tree lines off
hLines = findobj(hAxes,'Type','line');
set(hLines,'HandleVisibility','off');
set(hAxes,'TickLength',[0 0]);
set(hFig,'CurrentAxes',hAxes(1));
setappdata(hFig,'ClustergramStructure',clustStruct);

localSetupLimitListener(hAxes(1))
localSetupPositionListener(hAxes(1),hAxes(2:end));
localSetupImageDeleteListener(hAxes(1));


if nargout
    varargout{1} = hAxes;
    varargout{2} = hImage;
    varargout{3} = hDen;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, rowLabels, colLabels ] = makeRandomData
%makeRandomData creates some random data

data = sort(randn(26,4));
data = data(randperm(size(data,1)),:);

rowLabels = cellstr([ repmat('Row ',26,1) char(('A' + (0:25))')]);
colLabels = {'Col One','Col Two','Col Three','Col Four'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localMotionCallback(varargin)

hFig = localCleanUpLabel(varargin{1});

imageHandle = varargin{3};
if hittest(hFig) ~= imageHandle
    return
end

htext = varargin{4};

% get the labels from the image's UserData
ud = get(imageHandle,'UserData');
colLabels = ud.colLabels;
rowLabels = ud.rowLabels;

% Cdata is the actual values of the image
cdata = get(imageHandle,'Cdata');

% get the handle for the axes
imageAxes = get(imageHandle,'parent');

% get the position on the image.
cpAct = get(imageAxes,'CurrentPoint');

% determine where in axes was clicked
x = cpAct(1,1);
y = cpAct(1,2);

% calculate indices into image CData
xind = round(x);
yind = round(y);

try
    % get the value and labels of the current point
    val = num2str(cdata(yind,xind));
    if isempty(rowLabels{yind})
        rowLabel = '';
    else
        rowLabel = rowLabels(yind);
    end
    if isempty(colLabels{xind})
        colLabel = '';
    else
        colLabel = colLabels(xind);
    end

    % create a new text object -- start with it invisible
    textstr = {};
    if ~isempty(val)
        textstr = [textstr {val}];
    end
    if ~isempty(colLabel)
        textstr = [textstr {char(colLabel)}];
    end
    if ~isempty(rowLabel)
        textstr = [textstr {char(rowLabel)}];
    end

    if isempty(htext) || ~ishandle(htext)
        % give it an off white background, black text and grey border
        htext = text(x,y,textstr,...
            'Visible','off',...
            'BackgroundColor', [1 1 0.933333],...
            'Color', [0 0 0],...
            'EdgeColor', [0.8 0.8 0.8],...
            'Tag','ImageDataTip',...
            'interpreter','none');
    end



    % determine the offset in pixels
    set(htext,'position',[x y]);
    xy = [x y];
    %meanxy = [mean(xdata) mean(ydata)];
    meanxy = [mean(xlim(imageAxes)) mean(ylim(imageAxes))];
    setTextLocation(htext,xy,meanxy)


    % show the text
    set(htext, 'Visible','on')
    set(hFig,'WindowButtonUpFcn',{@localToggleOff,imageHandle,htext});
catch
    % quietly do nothing if something is out of range
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = localCleanUpLabel(varargin)
%localCleanUpLabel callback function to remove label from image

% get the handles to the figure, image and axis
if nargin == 0
    hFig = gcbf;
elseif ishandle(varargin{1})
    if strcmp('figure',get(varargin{1},'Type'))
        hFig = varargin{1};
    elseif strcmp('image',get(varargin{1},'Type'))
        hAxes = get(varargin{1},'Parent');
        hFig = get(hAxes,'Parent');
    end
end

% delete the old label if it exists
oldLabel = findobj(hFig,'Tag','ImageDataTip','Type','text');
if ~isempty(oldLabel)
    delete(oldLabel);
end

if nargout
    varargout{1} = hFig;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function local1DImageButtonCallback(varargin)
% localImageButtonCallback
% shows labels and values at selected point.

% clean up any old labels
hFig = localCleanUpLabel(varargin{1});

imageHandle = findobj(hFig,'type','image','Tag','ClusterGramHeatMap');

% get the labels from the image's UserData
ud = get(imageHandle,'UserData');
colLabels = ud.colLabels;
rowLabels = ud.rowLabels;

% Cdata is the actual values of the image
cdata = get(imageHandle,'Cdata');

% get the handle for the axes
imageAxes = get(imageHandle,'parent');

% get the position on the image.
cpAct = get(imageAxes,'CurrentPoint');

% determine where in axes was clicked
x = cpAct(1,1);
y = cpAct(1,2);

% calculate indices into image CData
xind = round(x);
yind = round(y);

try
    % get the value and labels of the current point
    val = num2str(cdata(yind,xind));
    if isempty(rowLabels{yind})
        rowLabel = '';
    else
        rowLabel = rowLabels(yind);
    end
    if isempty(colLabels{xind})
        colLabel = '';
    else
        colLabel = colLabels(xind);
    end

    % create a new text object -- start with it invisible
    textstr = {};
    if ~isempty(val)
        textstr = [textstr {val}];
    end
    if ~isempty(colLabel)
        textstr = [textstr {char(colLabel)}];
    end
    if ~isempty(rowLabel)
        textstr = [textstr {char(rowLabel)}];
    end

    htext = text(x,y,textstr,'visible','off');

    % give it an off white background, black text and grey border
    set(htext, 'BackgroundColor', [1 1 0.933333],...
        'Color', [0 0 0],...
        'EdgeColor', [0.8 0.8 0.8],...
        'Tag','ImageDataTip',...
        'interpreter','none');


    % determine the offset in pixels
    set(htext,'position',[x y]);
    xy = [x y];
    % meanxy = [mean(xdata) mean(ydata)];
    meanxy = [mean(xlim(imageAxes)) mean(ylim(imageAxes))];
    setTextLocation(htext,xy,meanxy)


    % show the text
    set(htext, 'Visible','on')
    set(hFig,'WindowButtonMotionFcn',{@localMotionCallback,imageHandle,htext});
    set(hFig,'WindowButtonUpFcn',{@localToggleOff,imageHandle,htext});
catch
    % quietly do nothing if something is out of range
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function local2DImageButtonCallback(varargin)
% localImageButtonCallback
% shows labels and values at selected point.

% clean up any old labels
hFig = localCleanUpLabel(varargin{1});

imageHandle = findobj(hFig,'type','image','tag','ClusterGramHeatMap');

% get the labels from the image's UserData
ud = get(imageHandle,'UserData');
colLabels = ud.colLabels;
rowLabels = ud.rowLabels;

% Cdata is the actual values of the image
cdata = get(imageHandle,'Cdata');

% get the handle for the axes
imageAxes = get(imageHandle,'parent');

% get the position on the image.
cpAct = get(imageAxes,'CurrentPoint');

% determine where in axes was clicked
x = cpAct(1,1);
y = cpAct(1,2);

% calculate indices into image CData
xind = round(x);
yind = floor(y + .5);

try
    % get the value and labels of the current point
    val = num2str(cdata(yind,xind));
    if isempty(rowLabels{yind})
        rowLabel = '';
    else
        rowLabel = rowLabels(yind);
    end
    if isempty(colLabels{xind})
        colLabel = '';
    else
        colLabel = colLabels(xind);
    end

    % create a new text object -- start with it invisible
    textstr = {};
    if ~isempty(val)
        textstr = [textstr {val}];
    end
    if ~isempty(colLabel)
        textstr = [textstr {char(colLabel)}];
    end
    if ~isempty(rowLabel)
        textstr = [textstr {char(rowLabel)}];
    end

    htext = text(x,y,textstr,'visible','off');

    % give it an off white background, black text and grey border
    set(htext, 'BackgroundColor', [1 1 0.933333],...
        'Color', [0 0 0],...
        'EdgeColor', [0.8 0.8 0.8],...
        'Tag','ImageDataTip',...
        'interpreter','none');

    % determine the offset in pixels
    set(htext,'position',[x y]);
    xy = [x y];
    % meanxy = [mean(xdata) mean(ydata)];
    meanxy = [mean(xlim(imageAxes)) mean(ylim(imageAxes))];
    setTextLocation(htext,xy,meanxy)

    % show the text
    set(htext, 'Visible','on')
    set(hFig,'WindowButtonMotionFcn',{@localMotionCallback,imageHandle,htext});
    set(hFig,'WindowButtonUpFcn',{@localToggleOff,imageHandle,htext});
catch
    % quietly do nothing if something is out of range
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hFig] = localToggleOff(varargin)
%LOCALTOGGLEOFF callback function to remove disable dragging of labels

hFig = varargin{1};

% clean up the labels
localCleanUpLabel;

imageHandle =varargin{3};

ud = get(imageHandle,'UserData');
if ud.onedimension
    set(imageHandle,'ButtonDownFcn',@local1DImageButtonCallback);
else
    set(imageHandle,'ButtonDownFcn',@local2DImageButtonCallback);
end

set(hFig,'WindowButtonMotionFcn','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = create1DClustergram(clustStruct,hFig)

if clustStruct.optLeafOrder
    % calculate pairwise distances and linkage
    dist = pdist(clustStruct.data,clustStruct.pdistargs{:});
    Z1 = linkage(dist,clustStruct.linkageargs);
    order = optimalleaforder(Z1,dist);
    clear('dist');
    if isempty(order) || (numel(unique(order)) ~= numel(order))
        warning('Bioinfo:clustergram:optLeafOrderFailure',...
            'The optimal leaf ordering algorithm failed to converge. Using default ordering.')
        [H1, T1, perm1] = dendrogram(Z1,0,clustStruct.dendroargs{:},'orientation','left'); %#ok
    else
        [H1, T1, perm1] = dendrogram(Z1,0,clustStruct.dendroargs{:},'orientation','left','r',order); %#ok
    end
else
    % calculate pairwise distances and linkage without storing the distance
    % matrix
    Z1 = linkage(clustStruct.data,clustStruct.linkageargs,clustStruct.pdistargs);
    [H1, T1, perm1] = dendrogram(Z1,0,clustStruct.dendroargs{:},'orientation','left'); %#ok
end

dAxes1 = get(H1(1),'Parent');

% make sure the dendrogram is in the proper figure
dAxes1Parent = get(dAxes1,'Parent');
if  clustStruct.reuse && dAxes1Parent ~= hFig
    set(dAxes1,'Parent',hFig)
    delete(dAxes1Parent);
end

set(hFig,'DoubleBuffer','on',...
    'Tag','ClustergramFigure')

% Move the dendrogram to the left edge of the axes, make the axes larger
fullWidth = clustStruct.fullwidth;
fullHeight = clustStruct.fullheight;

dWidth = fullWidth * clustStruct.ratio(1);
if clustStruct.onedimension
    dHeight = fullHeight;
else
    dHeight = fullHeight * (1 - clustStruct.ratio(2));
end

dPos = [(1 - fullWidth)/2, (1 - fullHeight)/2, dWidth, dHeight];

set(dAxes1,'Position',dPos,...
    'Visible','off',...
    'YTick',[],...
    'XTick',[],...
    'HandleVisibility','off');

% use a pretty colormap
set(hFig,'Colormap',clustStruct.colormap);

% we assume that the data is centred around zero
if clustStruct.zeroCenter
    maxval = min(max(abs(clustStruct.data(:))),clustStruct.saturation);
    minval = -maxval;
else
    maxval = min(max(clustStruct.data(:)),clustStruct.saturation);
    minval = min(clustStruct.data(:));
end
% determine the size of the image
nrows = size(clustStruct.data,1);
ncols = size(clustStruct.data,2);

imWidth = fullWidth-dWidth;
imHeight = dHeight;
imPos = [dPos(1) + dPos(3),dPos(2),imWidth,imHeight];

imAxes = axes('Parent',hFig,...
    'Units','normalized',...
    'Position',imPos);

hImage1 = imagesc(clustStruct.data(perm1,:),...
    'Parent',imAxes,...
    'Tag','ClustergramImage',...
    [minval,maxval]);

% Visible set to on by IMAGESC. Set to 'off'
set(imAxes,'Visible','off',...
    'YDir','normal');

% make sure labels are in a column cell array
if size(clustStruct.rowLabels,2) ~= 1
    clustStruct.rowLabels = clustStruct.rowLabels';
end

% store the labels in the UserData of the image
ud.rowLabels = clustStruct.rowLabels(perm1);
ud.colLabels = clustStruct.colLabels;

ud.onedimension = clustStruct.onedimension;

set(hImage1,'UserData',ud);

% put labels on the rows
if nrows < 200
    ytickvals = 1:nrows;
else
    ytickvals = [];
end

% put labels on the columns
if ncols < 200
    xtickvals =  1:ncols;
else
    xtickvals = [];
end

% get rid of the original row labels
t = findall(dAxes1,'Type','text','Tag','ClustergramRowLabel');
if ~isempty(t)
    delete(t)
end

% create text objects for labels
set(imAxes,'YTick',ytickvals,...
    'YTickLabel',clustStruct.rowLabels(perm1),...
    'XTick',xtickvals,...
    'XTickLabel',clustStruct.colLabels,...
    'Visible','on',...
    'YAxisLocation','right');

setFontSize(imAxes);

% now wouldn't it be nice if we could click on the image to and get some
% information about the particular spot.

set(hImage1,'ButtonDownFcn',@local1DImageButtonCallback,'tag','ClusterGramHeatMap')

set(dAxes1,'NextPlot','replace')

bh = hggetbehavior(dAxes1,'Zoom');
set(bh,'Enable',false);

hYLimLink = linkprop([imAxes dAxes1],'YLim');
setappdata(imAxes,'ClustergramYLimLink',hYLimLink)

if nargout
    varargout{1} = [imAxes dAxes1];
    varargout{2} = hImage1;
    varargout{3} = H1;
    varargout{4} = perm1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = create2DClustergram(clustStruct,hFig1)

% create 1D first, to get axes, dendrogram and permutation index for data
[dAxes1,hImage1,H1,perm1] = create1DClustergram(clustStruct,hFig1);

imAxes = dAxes1(1);
dAxes1 = dAxes1(2);

% create a new, invisible figure to hold the second dendrogram, before moving it to the
% existing figure
cgfig2 = figure('Visible','off');

if clustStruct.optLeafOrder
    % calculate pairwise distances and linkage
    dist = pdist(clustStruct.data(perm1,:)',clustStruct.pdistargs{:});
    Z2 = linkage(dist,clustStruct.linkageargs);
    order = optimalleaforder(Z2, dist);
    clear('dist');
    if isempty(order) || (numel(unique(order)) ~= numel(order))
        warning('Bioinfo:clustergram:optLeafOrderFailure',...
            'The optimal leaf ordering algorithm failed to converge. Using default ordering.')
        [H2, T2, perm2] = dendrogram(Z2,0,clustStruct.dendroargs{:},'orientation','top'); %#ok
    else

        [H2, T2, perm2] = dendrogram(Z2,0,clustStruct.dendroargs{:},'orientation','top','r',order); %#ok
    end
else
    % calculate pairwise distances and linkage without storing the distance
    % matrix
    Z2 = linkage(clustStruct.data(perm1,:)',clustStruct.linkageargs,clustStruct.pdistargs);
    [H2, T2, perm2] = dendrogram(Z2,0,clustStruct.dendroargs{:},'orientation','top'); %#ok
end
% % calculate pairwise distances and linkage
% Z2 = linkage(clustStruct.data(perm1,:)',clustStruct.linkageargs,clustStruct.pdistargs);
%
%
% [H2, T2, perm2] = dendrogram(Z2,0,clustStruct.dendroargs{:},'orientation','top'); %#ok
hiddenAxes = get(H2(1),'Parent');

hFig1 = get(dAxes1(1),'Parent');
imPos = get(imAxes,'Position');

fullHeight = clustStruct.fullheight;
dHeight = fullHeight * clustStruct.ratio(2);
dPos = [imPos(1),imPos(2) + imPos(4),imPos(3), dHeight];
dAxes2 = axes('Parent',hFig1,...
    'Units','normalized',...
    'Position',dPos,...
    'Visible','off',...
    'XLim',get(hiddenAxes,'XLim'),...
    'XTick',[],...
    'YTick',[],...
    'HandleVisibility','off');



% move the dendrogram lines into the original axes
set(H2,'Parent',dAxes2)

% delete the invisible figure
delete(cgfig2);

% delete the existing image
delete(hImage1);

% use a pretty colormap
colormap(clustStruct.colormap);

% we assume that the data is centred around zero
if clustStruct.zeroCenter
    maxval = min(max(abs(clustStruct.data(:))),clustStruct.saturation);
    minval = -maxval;
else
    maxval = min(max(clustStruct.data(:)),clustStruct.saturation);
    minval = min(clustStruct.data(:));
end

% determine the size of the image
nrows = size(clustStruct.data,1);
ncols = size(clustStruct.data,2);

hImage2 = imagesc(clustStruct.data(perm1,perm2),...
    'Parent',imAxes,...
    'Tag','ClustergramImage',...
    [minval,maxval]);

set(imAxes,'YDir','normal','Visible','off');

% make sure labels are in a column cell array
if size(clustStruct.rowLabels,2) ~= 1
    clustStruct.rowLabels = clustStruct.rowLabels';
end

% store the labels in the UserData of the image
ud.rowLabels = clustStruct.rowLabels(perm1);
ud.colLabels = clustStruct.colLabels(perm2);

ud.onedimension = clustStruct.onedimension;

set(hImage2,'UserData',ud);

% fixDendrograms(H1,H2,hImage2,clustStruct.ratio);

% put labels on the rows
if nrows < 200
    ytickvals = 1:nrows;
else
    ytickvals = [];
end


% put labels on the columns
if ncols < 200
    xtickvals =  1:ncols;
else
    xtickvals = [];
end

% create text objects for labels
set(imAxes,'YTick',ytickvals,...
    'YTickLabel',clustStruct.rowLabels(perm1),...
    'XTick',xtickvals,...
    'XTickLabel',clustStruct.colLabels(perm2),...
    'Visible','on',...
    'YAxisLocation','right');

setFontSize(imAxes)
% now wouldn't it be nice if we could click on the image to and get some
% information about the particular spot.

set(hImage2,'ButtonDownFcn',@local2DImageButtonCallback,'tag','ClusterGramHeatMap')

set(dAxes1,'NextPlot','replace')

bh = hggetbehavior(dAxes2,'Zoom');
set(bh,'Enable',false);

hXLimLink = linkprop([imAxes dAxes2],'XLim');
setappdata(imAxes,'ClustergramXLimLink',hXLimLink);

if nargout
    varargout{1} = [imAxes dAxes1 dAxes2];
    varargout{2} = hImage2;
    varargout{3} = [H1;H2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setTextLocation(htext,xy,meanxy)

textUnits = get(htext,'Units');
set(htext,'Units','pixels')

pixpos = get(htext,'Position');

offsets = [0 0 0];

% determine what quadrant the pointer is in
quadrant= xy < meanxy;

if ~quadrant(1)
    set(htext,'HorizontalAlignment','Right')
    offsets(1) = -4;
else
    set(htext,'HorizontalAlignment','Left')
    offsets(1) = 16;
end
if quadrant(2)
    set(htext,'VerticalAlignment','Bottom')
    offsets(2) = 4;
else
    set(htext,'VerticalAlignment','Top')
    offsets(2) = -4;
end

set(htext,'Position',pixpos + offsets);
set(htext,'Units',textUnits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localSetupLimitListener(imAxes)
% helper function to setsup listeners for the ylables, so we can detect if
% we would need to change the fontsize
hgp     = findpackage('hg');
axesC   = findclass(hgp,'axes');
LimListener = handle.listener(imAxes,[axesC.findprop('XLim') axesC.findprop('YLim')],...
    'PropertyPostSet',@localLimitListener);

hFig = get(imAxes,'Parent');
listeners = getappdata(hFig,'ClustergramListeners');

setappdata(hFig,'ClustergramListeners',[listeners LimListener]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localLimitListener(hSrc,event)%#ok

imAxes = event.AffectedObject;
newLim = get(event,'NewValue');
hImage = findobj(imAxes,'Type','image');
ud = get(hImage,'UserData');

if strcmp(get(event.Source,'Name'),'XLim') %XLim
    if diff(newLim) < 200 && ~isempty(ud) && isstruct(ud)
        xtickvals = max(1,ceil(newLim(1))):min(floor(newLim(2)),numel(ud.colLabels));
        xticklabels =  ud.colLabels(xtickvals);
    else
        xtickvals = [];
        xticklabels = '';
    end
    set(imAxes,'XTick',xtickvals,'XTickLabel',xticklabels);
else %YLim
    if diff(newLim) < 200 && ~isempty(ud) && isstruct(ud)
        ytickvals = max(1,ceil(newLim(1))):min(floor(newLim(2)),numel(ud.rowLabels));
        yticklabels =  ud.rowLabels(ytickvals);
    else
        ytickvals = [];
        yticklabels = '';
    end
    set(imAxes,'YTick',ytickvals,'YTickLabel',yticklabels);
end
setFontSize(imAxes);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localSetupPositionListener(imAxes,dendroAxes)

hgp     = findpackage('hg');
axesC   = findclass(hgp,'axes');
PostPositionListener = handle.listener(imAxes,axesC.findprop('Position'),...
    'PropertyPostSet',{@localPostPositionListener,dendroAxes});

hFig = get(imAxes,'Parent');
listeners = getappdata(hFig,'ClustergramListeners');

setappdata(hFig,'ClustergramListeners',[listeners PostPositionListener]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localPostPositionListener(hSrc,event,dAxes) %#ok

imAxes = event.AffectedObject;
fixAxesLocations(imAxes,dAxes);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localSetupImageDeleteListener(imAxes)

ImageDeleteListener = handle.listener(imAxes,'ObjectChildRemoved',@localImageDeleteListener);

hFig = get(imAxes,'Parent');
listeners = getappdata(hFig,'ClustergramListeners');

setappdata(hFig,'ClustergramListeners',[listeners ImageDeleteListener]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function localImageDeleteListener(hSrc,event)

if isa(event.Child,'image')
    hFig = get(hSrc,'Parent');
    clf(hFig,'reset')
    rmappdata(hFig,'ClustergramListeners');
    rmappdata(hFig,'ClustergramStructure');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fixAxesLocations(imAxes,dAxes)


imUnits = get(imAxes,'Units');
dAxes1Units = get(dAxes(1),'Units');

set(imAxes,'Units','pixels');
set(dAxes(1),'Units','pixels');

imPos = get(imAxes,'Position');
dAxes1Pos = get(dAxes(1),'Position');

% match edges
if dAxes1Pos(1) + dAxes1Pos(3) ~= imPos(1)
    dAxes1Pos(1) = imPos(1) - dAxes1Pos(3);
end

% y
if dAxes1Pos(2) ~= imPos(2)
    dAxes1Pos(2) = imPos(2);
end

% height
if dAxes1Pos(4) ~= imPos(4)
    dAxes1Pos(4) = imPos(4);
end

set(dAxes(1),'Position',dAxes1Pos);
set(dAxes(1),'Units',dAxes1Units);

% second axes
if numel(dAxes) == 2

    dAxes2Units = get(dAxes(2),'Units');
    set(dAxes(2),'Units','pixels');
    dAxes2Pos = get(dAxes(2),'Position');

    % match edges
    if dAxes2Pos(2) ~= imPos(2) + imPos(4)
        dAxes2Pos(2) = imPos(2) + imPos(4);
    end

    % x
    if dAxes2Pos(1) ~= imPos(1)
        dAxes2Pos(1) = imPos(1);
    end

    % width
    if dAxes2Pos(3) ~= imPos(3)
        dAxes2Pos(3) = imPos(3);
    end

    set(dAxes(2),'Position',dAxes2Pos);
    set(dAxes(2),'Units',dAxes2Units);

end

set(imAxes,'Units',imUnits)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setFontSize(imAxes)
% Try to keep font size reasonable for x and y tick labels
hFig = get(imAxes,'Parent');
magicNumber = 200;
nrows = diff(get(imAxes,'YLim'));
ncols = diff(get(imAxes,'XLim'));
if ncols < magicNumber && nrows < magicNumber
    ratio = max(get(hFig,'Position').*[0 0 0 1])/max(nrows,ncols);
elseif ncols < magicNumber
    ratio = max(get(hFig,'Position').*[0 0 0 1])/ncols;
elseif nrows < magicNumber
    ratio = max(get(hFig,'Position').*[0 0 0 1])/nrows;
else
    ratio = 1;
end
set(imAxes,'Fontsize',min(9,ceil(ratio/1.5)));    % the gold formula
