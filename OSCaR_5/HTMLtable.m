function HTML = HTMLtable(matrix,varargin)
% HTMLtable - Generate the code for an CSS/HTML-based data table from an input numeric, string, or cell array.
% The generated HTML code is compatible with app designer as the 'HTMLSource' for an HTML element.
% With 60 user adjustable parameters, the table can be customized in
% virtually any way.
%
% Author: Eric Dziekonski, Ph.D.
% Email: etdzi1@gmail.com
% Date: March 30, 2020
% Version: 1.0.1
%
% Description of arguments:
% Note: All name-value pairs are optional.
%
% matrix (REQUIRED): Value must be a numeric or cell array. Cells must 
%       contain numeric, character, or string values. Nested tables are not
%       supported.
% 'PageBackgroundColor': Background color of the webpage. Valid value:
%       [1x3] rgb color scaled between 0 and 255.
% 'Title': Table title. Valid value: 1D character array or [1x1] string array
% 'TitleFontName': Font name of the table title. Valid values: any available 
%       system font (see listfonts).
% 'TitleFontSize': Font size of table title in pixels. Valid value: positive 
%       scalar integer.
% 'TitleFontWeight': Font weight of table title. Valid values: 'normal',
%       'bold','bolder','lighter'.
% 'TitleFontStyle': Style of table title. Valid values: 'normal','italic',
%       'oblique'.
% 'TitleFontColor': Title text color. Valid value: [1x3] rgb color scaled 
%       between 0 and 255.
% 'TitleTextAlign': Alignment of table title. Valid values: 'left','center',
%       'right'.
% 'TitleBackgroundColor': Background color of the table title area. Valid 
%       value: [1x3] rgb color scaled between 0 and 255.
% 'LabelLayer': Include an additional header for tables that have more than 
%       two dimensions. The header indicates the layer of the table being
%       displayed. If, for instance, you have a table that is [1000x20] and
%       you want to display column headings every 100 rows, you can reshape
%       the table to [100x20x10] and make 'LabelLayer' false as they are
%       all originally of a single layer. Valid values: 'true','false'.
% 'ColumnNames': [1xm] or [mx1] numeric or string array, or a [1xm] or [mx1] 
%       cell array of numeric, character, or string values for the column
%       headings, where 'm' is the number of columns in the matrix. All 
%       parameters starting with 'Column' act on the value of 'ColumnNames'.
%       Therefore, if 'ColumnNames' is not populated, parameters starting
%       with 'Column' will have no effect on the generated table.
% 'ColumnFontName': Font name of all column headings. Valid values: any 
%       available system font (see listfonts).
% 'ColumnFontSize': Font size of column headings in pixels. Valid value: 
%       positive scalar integer.
% 'ColumnPadding': Padding of all column headings in pixels. Valid value: 
%       positive scalar integer.
% 'ColumnFormat': Value must be a valid sprintf output format, e.g. '%.3f'. 
%       Valid input formats: 1D char array (applied to all column headings)
%       [1x1] cell array (applied to all column headings), [1xm] or [mx1] 
%       cell array, where ''m'' is the number of columns in the data matrix 
%       (individual column heading formats).
% 'ColumnFontWeight': Font weight of column headings. Valid values: 'normal',
%       'bold','bolder','lighter'.
% 'ColumnFontStyle': Style of column headings. Valid values: 'normal','italic',
%       'oblique'.
% 'ColumnFontColor': Value must be a [1x3] numeric array (applied to all 
%       column headings) or a [mx3] numeric array (applied to column headings 
%       individually), where 'm' is the number of columns in the data matrix. 
%       Values must be scaled between 0 and 255.
% 'ColumnMinimumWidth': Minimum column width in pixels. Valid values: 
%       positive scalar integer.
% 'ColumnTextAlign': Alignment of column headings. Valid values: 'left',
%       'center','right'.
% 'ColumnBackgroundColor': Value must be a [1x3] numeric array (applied to 
%       all column headings) or a [mx3] numeric array (applied to column 
%       headings individually), where 'm' is the number of columns in the 
%       data matrix. Values must be scaled between 0 and 255.
% 'RowNames': [1xn] or [nx1] numeric or string array, or a [1xn] or [nx1] cell array
%       of numeric, character, or string values for the row headings, where 'n' 
%       is the number of rows in the matrix. All parameters starting with 'Row' 
%       act on the value of 'RowNames'. Therefore, if 'RowNames' is not populated, 
%       parameters starting with 'Row' will have no effect on the generated 
%       table.
% 'RowFontName': Font name of all row headings. Valid values: any available 
%       system font (see listfonts).
% 'RowFontSize': Font size of row headings in pixels. Valid value: positive 
%       scalar integer.
% 'RowPadding': Padding of all row headings in pixels. Valid value: positive 
%       scalar integer.
% 'RowFormat': Value must be a valid sprintf output format, e.g. '%.3f'. 
%       Valid input formats: 1D char array (applied to all row headings)
%       [1x1] cell array (applied to all row headings), [1xm] or [mx1] 
%       cell array, where 'm' is the number of columns in the data matrix 
%       (individual row heading formats).
% 'RowFontWeight': Font weight of row headings. Valid values: 'normal',
%       'bold','bolder','lighter'.
% 'RowFontStyle': Style of row headings. Valid values: 'normal','italic',
%       'oblique'
% 'RowFontColor': Value must be a [1x3] numeric array (applied to all row 
%       headings) or a [nx3] numeric array (applied to row headings individually), 
%       where 'n' is the number of rows in the data matrix. Values must be 
%       scaled between 0 and 255.
% 'RowTextAlign': Alignment of row headings. Valid values: 'left','center',
%       'right'
% 'RowBackgroundColor': Value must be a [1x3] numeric array (applied to all 
%       row headings) or a [nx3] numeric array (applied to row headings 
%       individually), where ''n'' is the number of rows in the data matrix. 
%       Values must be scaled between 0 and 255.
% 'DataFontName': Font name of all table cells. Valid values: any available 
%       system font (see listfonts).
% 'DataFontSize': Font size of table cells in pixels. Valid value: positive 
%       scalar integer.
% 'DataPadding': Padding of all table cells in pixels. Valid value: positive 
%       scalar integer.
% 'DataFormat': Value must be a valid sprintf output format, e.g. ''%.3f''. 
%       Valid input formats: 1D char array (applied to all data), [1x1] cell 
%       array (applied to all data), [1xm] or [mx1] cell array, where ''m'' 
%       is the number of columns in the data matrix (column-dependent data 
%       format).
% 'DataFontWeight': Font weight of table cells. Valid values: 'normal','bold',
%       'bolder','lighter'.
% 'DataFontStyle': Style of table cells. Valid values: 'normal','italic',
%       'oblique'.
% 'DataFontColor': Value must be a [1x3], [nx3], [mx3], [nxmx3], or [n x m x 3k]
%       numeric array, where ''n'' is the number of rows in the data matrix, 
%       ''m'' is the number of columns in the data matrix, and ''k'' is the 
%       number of layers in the data matrix. Values must be scaled between 
%       0 and 255. In the case of n = m, set the preferred color direction 
%       ('DataFontColorDimension') to 'row' or 'column' (default). 
%       Alternatively, if the matrix is numeric, or a cell array of numeric 
%       values, one can specify a supported colormap (e.g. 'parula') to 
%       scale the values to. The data will be scaled along the dimension 
%       specified by 'DataFontColorDimension'.
% 'DataFontColorDimension': When 'DataFontColor' is a [1x3], [nx3], [mx3], 
%       [nxmx3], or [n x m x 3k], where n ~= m, this parameter does
%       nothing. If however, n == m and one specifies an rgb array of size 
%       [nx3] == [mx3] it becomes ambiguous and unclear if the user wants
%       the specified colors to be applied to the column or to the row.
%       Therefore, the user can specify the direction. This also holds true
%       if the user specifies a colormap. One can choose if the colormap is
%       to be scaled along the 'row' (n), 'column' (m), 'layer' (k), or 
%       the entire 'table'.
% 'NegativeDataFontColor': Font color of negative cell values. Reuires the 
%       input matrix to be numeric. Valid value: [1x3] rgb color scaled 
%       between 0 and 255.
% 'DataTextAlign': Alignment of table cells. Valid values: 'left','center',
%       'right'.
% 'DataBackgroundColor': Value must be a [1x3], [nx3], [mx3], [nxmx3], or [n x m x 3k]
%       numeric array, where ''n'' is the number of rows in the data matrix, 
%       ''m'' is the number of columns in the data matrix, and ''k'' is the 
%       number of layers in the data matrix. Values must be scaled between 
%       0 and 255. In the case of n = m, set the preferred color direction 
%       ('DataBackgroundColorDimension') to 'row' or ''column' (default). 
%       Alternatively, if the matrix is numeric, or a cell array of numeric 
%       values, one can specify a supported colormap (e.g. 'parula') to 
%       scale the values to. The data will be scaled along the dimension 
%       specified by 'DataBackgroundColorDimension'. While 'RowStriping' is 
%       enabled, 'DataBackgroundColor' will NOT be applied to the table cells.
% 'DataBackgroundColorDimension': When 'DataBackgroundColor' is a [1x3], [nx3], 
%       [mx3], [nxmx3], or [n x m x 3k], where n ~= m, this parameter does
%       nothing. If however, n == m and one specifies an rgb array of size 
%       [nx3] == [mx3] it becomes ambiguous and unclear if the user wants
%       the specified colors to be applied to the column or to the row.
%       Therefore, the user can specify the direction. This also holds true
%       if the user specifies a colormap. One can choose if the colormap is
%       to be scaled along the 'row' (n), 'column' (m), 'layer' (k), or 
%       the entire 'table'.
% 'DataEditable': When enabled, table cells are editable. The new value can 
%       be retrieved by enabling 'RecordRowClick'. Valid values: 'true','false'.
% 'OmitValues': When false, cell text is not populated or displayed in the table. 
%       Therefore, the cells themselves are not editable. Another approach
%       is to have the data cell background color and the cell font color
%       to be the same. This hides the data from view, but the table is
%       still populated. In this approach, the table can still be edited
%       and the values can still be viewed if 'HighlightHover' is active.
%       Valid values: 'true','false'.
% 'Striping': Enable/disable row striping. While enabled, 'DataBackgroundColor'
%       will NOT be applied to the table cells. Valid values: 'true','false'.
% 'StripingColors': Value must be a [2x3] numeric array with values scaled 
%       between 0 and 255. The color defined in row one is applied to odd 
%       rows of the table. The color defined in row two is applied to even 
%       rows of the table.
% 'BorderStyle': Style of the border. To turn the border off, set the value 
%       to 'none'. If set to 'none', parameters starting with 'Border' will 
%       have no effect on the table. Valid values: 'none','hidden','dotted',
%       'dashed','solid','double','groove','ridge','inset','outset'. The
%       styles of 'groove', 'ridge', 'inset', and 'outset' may require you
%       to change the border color/width in order to visualize the effect.
% 'BorderWidth': Width of the border in pixels. Valid value: positive scalar 
%       integer.
% 'BorderColor': Border color. Valid value: [1x3] rgb color scaled between 
%       0 and 255.
% 'BorderCollapse': Sets whether table borders should collapse into a single 
%       border or be separated as in standard HTML. Valid values: 'separate',
%       'collapse'.
% 'HighlightHover': Highlight the current row/cell when hovered over. 
%       Valid values: 'true','false'.
% 'HoverColor': Color of cell(s) when the mouse is hovering over them. Valid 
%       value: [1x3] rgb color scaled between 0 and 255.
% 'HighlightRows': Value must be a 1D numeric array containing the row numbers 
%       to be highlighted. Negative numbers are ignored. Values are rounded 
%       to the nearest integer.
% 'HighlightColor': Value must be a [1x3] numeric array or [nx3] numeric array, 
%       where ''n'' is the number of values in the numeric array of 
%       ''HighlightRows''. Values must be scaled between 0 and 255.
% 'ShowOutput': When true, the function opens a webpage to view the resulting 
%       HTML table. Valid values: 'true','false'.
% 'RecordRowClick': If enabled, a three value numeric array is written to the 
%       objects 'Data' property. It contains:
%       [row number clicked on, last row clicked on, value in last cell clicked on].
%       One can use this to implement a 'DataChangedFcn' for the table when
%       used as the HTMLSource for a uihtml element in app designer, for example.
%       The 'value in the last cell clicked on' is the edited value, if
%       applicable. Therefore, one can use the 'DataChangedFcn' to check if
%       the cell data was edited. If so, one could recall the HTMLtable
%       function to update cell colors accordingly using the new cell
%       value. Valid values: 'true','false'.
% 'FileName': Output file name. '.html' is not necessary. The path may be
%       included to save to another directory.
% 'OutputToFile': When true, the function saves the generated HTML text to 
%       the file specified by 'FileName'. A valid 'FileName' is required. 
%       Valid values: 'true','false'.
%
% Tips/Tricks:
% 1) You can modify the default table behavior by adjusting the default
%       structure below. You may do this to achieve the desired look, after 
%       which the function can be called with minimal name-value pairs. If 
%       you do this, pay careful attention to the value formats and their 
%       dimensions as this bypasses the input parser which checks the input 
%       format.
% 2) The 'name' parameter is case insensitive and allows for partial matching.
%       The 'value' parameter of the name-value pairs is case insensitive but 
%       requires the full value. This was done to avoid using
%       'validatestring' in the parser which took a significant amount of
%       time. Instead, 'strcmpi' checks the validity of the value.
% 3) With small tables, a significant portion of the functions execution time
%       is due to 'listfonts'. One may copy Matlabs built-in listfonts function,
%       ('open listfonts'), and include it here as a subfunction. Comment out 
%       the final two lines of the function which sort the list and removes 
%       duplicates; this feature is uneccessary and can shave off 5-7ms.
% 4) Faster execution times can be achieve with:
%       'DataEditable',false
%       'ShowOutput',false
%       'OutputToFile',false
%       'RowStriping',true
%
% Acknowledgements:
% Inspiration for additional features was drawn from:
% Gus Brown (2020). Display data as an HTML table 
%       https://www.mathworks.com/matlabcentral/fileexchange/18329-display-data-as-an-html-table, 
%       MATLAB Central File Exchange. Retrieved March 24, 2020.
% Roger Parkyn (2020). Html Table Writer 
%       https://www.mathworks.com/matlabcentral/fileexchange/25078-html-table-writer, 
%       MATLAB Central File Exchange. Retrieved March 24, 2020.
%
% % EXAMPLES:
% 
% % Example input matricies
% Matrix1 = rand(10,5,2); % 3D numerical
% Matrix1(:,:,1) = Matrix1(:,:,1) .* 0.5; % Change scaling of layer 1
% Matrix2 = string(Matrix1); % 3D string array
% Matrix3 = num2cell(Matrix1); % 3D cell array of numerical values
% Matrix4 = cellfun(@num2str,num2cell(Matrix1),'UniformOutput',false); % 3D cell array of character values
% Matrix5 = cellfun(@string,cellfun(@num2str,num2cell(Matrix1),'UniformOutput',false),'UniformOutput',false); % 3D cell array of string values
%
% Title = 'Example Table';
% RowNames = (1:size(Matrix1,1)).';
% ColumnNames = cellfun(@num2str,num2cell(1:size(Matrix1,2)),'UniformOutput',false);
% 
% % A basic example:
% out = HTMLtable(Matrix1(:,:,1),'ShowOutput',true);  
%
% % Example with all properties:
% htmlcode = HTMLtable(...
%     Matrix1,... 
%     'PageBackgroundColor',[],...
%     'Title',Title,...                          
%     'TitleFontName','Times New Roman',...   
%     'TitleFontSize',20,...                  
%     'TitleFontWeight','bold',...           
%     'TitleFontStyle','italic',...          
%     'TitleFontColor',[],...                
%     'TitleTextAlign','center',...         
%     'TitleBackgroundColor',[],...         
%     'LabelLayer',true,...
%     'ColumnNames',ColumnNames,...             
%     'ColumnFontName','Helvetica',...      
%     'ColumnFontSize',9,...                 
%     'ColumnPadding',3,...                  
%     'ColumnFormat','',...                  
%     'ColumnFontWeight','bold',...        
%     'ColumnFontStyle','normal',...          
%     'ColumnFontColor',[],...               
%     'ColumnMinimumWidth',100,...            
%     'ColumnTextAlign','center',...          
%     'ColumnBackgroundColor',[],...          
%     'RowNames',RowNames,...                 
%     'RowFontName','Helvetica',...         
%     'RowFontSize',9,...                    
%     'RowPadding',3,...                      
%     'RowFormat','',...                      
%     'RowFontWeight','bold',...            
%     'RowFontStyle','normal',...             
%     'RowFontColor',[],...                  
%     'RowTextAlign','center',...             
%     'RowBackgroundColor',[],...            
%     'DataFontName','Helvetica',...          
%     'DataFontSize',9,...                    
%     'DataPadding',3,...                     
%     'DataFormat','',...                     
%     'DataFontWeight','normal',...           
%     'DataFontStyle','normal',...            
%     'DataFontColor',[],...                 
%     'DataFontColorDimension','',...      
%     'NegativeDataFontColor',[255,100,100],...          
%     'DataTextAlign','left',...              
%     'DataBackgroundColor',[],...            
%     'DataBackgroundColorDimension','',...
%     'DataEditable',true,...
%     'OmitValues',false,...
%     'Striping',true,...                   
%     'StripingColors',[],...               
%     'BorderStyle','solid',...              
%     'BorderWidth',1,...                     
%     'BorderColor',[],...                   
%     'BorderCollapse','separate',...         
%     'HighlightHover',true,...              
%     'HoverColor',[],...                    
%     'HighlightRows',[1,5,10],...                 
%     'HighlightColor',[100,255,100],...                 
%     'ShowOutput',true,...                  
%     'RecordRowClick',true,...               
%     'FileName','testFile',...                      
%     'OutputToFile',false);     
%
%
% % APP DESIGNER EXAMPLE
% % Set 'RecordRowClick' to true to enable the 'DataChangedFcn'
% fig = uifigure;
% h1 = uihtml(fig,...
%     'Position',[0,0,fig.Position(3),fig.Position(4)],...
%     'DataChangedFcn','output = h1.Data;');
% drawnow;
%
% h1.HTMLSource = htmlcode;

HTML = '';

%======================================================
% Define default values
%======================================================

default.pagebackgroundcolor = [255,255,255];
default.title = '';
default.titlefontname = 'Helvetica';
default.titlefontsize = 20;
default.titlefontweight = 'bold';
default.titlefontstyle = 'normal';
default.titlefontcolor = [0,0,0];
default.titletextalign = 'center';
default.titlebackgroundcolor = [255,255,255];
default.labellayer = true;
default.columnnames = {};
default.columnfontname = 'Helvetica';
default.columnfontsize = 9;
default.columnpadding = 3;
default.columnformat = '';
default.columnfontweight = 'bold';
default.columnfontstyle = 'normal';
default.columnfontcolor = [0,0,0];
default.columnminimumwidth = 100;
default.columntextalign = 'center';
default.columnbackgroundcolor = [210,210,210];
default.rownames = {};
default.rowfontname = 'Helvetica';
default.rowfontsize = 9;
default.rowpadding = 3;
default.rowformat = '';
default.rowfontweight = 'bold';
default.rowfontstyle = 'normal';
default.rowfontcolor = [0,0,0];
default.rowtextalign = 'center';
default.rowbackgroundcolor = [210,210,210];
default.datafontname = 'Helvetica';
default.datafontsize = 9;
default.datapadding = 3;
default.dataformat = '';
default.datafontweight = 'normal';
default.datafontstyle = 'normal';
default.datatextalign = 'center';
default.databackgroundcolor = [255,255,255];
default.databackgroundcolordimension = 'column';
default.datafontcolor = [0,0,0];
default.datafontcolordimension = 'column';
default.negativedatafontcolor = [];
default.dataeditable = false;
default.omitvalues = false;
default.borderstyle = 'solid';
default.borderwidth = 1;
default.bordercolor = [0,0,0];
default.bordercollapse = 'collapse';
default.striping = true;
default.stripingcolors = [255,255,255;230,230,230];
default.highlighthover = true;
default.hovercolor = [153,217,255];
default.highlightrows = [];
default.highlightcolor = [100,255,100];
default.showoutput = false;
default.recordrowclick = true;
default.filename = '';
default.outputtofile = false;

%======================================================
% Define valid expected values
%======================================================

expectedWeights = {'normal','bold','bolder','lighter'};
expectedStyles = {'normal','italic','oblique'};
expectedAlignments = {'left','center','right'};
expectedFonts = listfonts;
expectedBorderStyles = {'none','hidden','dotted','dashed','solid','double','groove','ridge','inset','outset'};
expectedBorderCollapse = {'separate','collapse'};
expectedDimensions = {'row','column','layer','table'};
expectedColormaps = {'parula','jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink','cividis','inferno','magma','plasma','twilight','viridis'};

%======================================================
% Parse inputs
%======================================================

p = inputParser;
p.CaseSensitive = false;
p.KeepUnmatched = true;
p.PartialMatching = true;
p.StructExpand = true;

checkScalar = @(x) assert(isempty(x) || (isnumeric(x) && isscalar(x) && (x > 0)),'Value must be positive, numeric, and scalar. An empty input sets the default value.');
checkLogicalScalar = @(x) assert(isempty(x) || (isscalar(x) && islogical(x)),sprintf('Expected a scalar logical value. The value was %s instead.',class(x)));
checkLogical1D = @(x) assert(isempty(x) || (islogical(x) && (length(x) == 1 || all(size(x,1,2) == [1,size(matrix,2)],'all') || all(size(x,1,2) == [1,size(matrix,2)],'all'))),sprintf('Value must be a [1x1] logical value or a [1xm] or [mx1] logical array, where ''m'' is the number of columns in the data matrix. The value was %s instead.',class(x)));
checkColor1D = @(x) assert(isempty(x) || (isnumeric(x) && all(size(x)==[1,3],'all') && all(x<=255 & x>=0,'all')),'Value must be a [1x3] numeric array of RGB values scaled between 0 and 255.');
checkString1D = @(x) assert(isempty(x) || ((ischar(x) || (isstring(x) && length(x) == 1)) && length(x)==numel(x)),'Value must be a 1-D character or [1x1] string array.');
checkFonts = @(x) assert(isempty(x) || any(strcmpi(x,expectedFonts)),sprintf(['Value expected to be one of the following: ' repmat('%s, ',1,length(expectedFonts)-1) '%s. The value was ''%s'' instead.'],expectedFonts{:},x)); 
checkWeights = @(x) assert(isempty(x) || any(strcmpi(x,expectedWeights)),sprintf(['Value expected to be one of the following: ' repmat('%s, ',1,length(expectedWeights)-1) '%s. The value was ''%s'' instead.'],expectedWeights{:},x)); 
checkStyle = @(x) assert(isempty(x) || any(strcmpi(x,expectedStyles)),sprintf(['Value expected to be one of the following: ' repmat('%s, ',1,length(expectedStyles)-1) '%s. The value was ''%s'' instead.'],expectedStyles{:},x)); 
checkAlignments = @(x) assert(isempty(x) || any(strcmpi(x,expectedAlignments)),sprintf(['Value expected to be one of the following: ' repmat('%s, ',1,length(expectedAlignments)-1) '%s. The value was ''%s'' instead.'],expectedAlignments{:},x)); 
checkBorderStyles = @(x) assert(isempty(x) || any(strcmpi(x,expectedBorderStyles)),sprintf(['Value expected to be one of the following: ' repmat('%s, ',1,length(expectedBorderStyles)-1) '%s. The value was ''%s'' instead.'],expectedBorderStyles{:},x)); 
checkBorderCollapse = @(x) assert(isempty(x) || any(strcmpi(x,expectedBorderCollapse)),sprintf(['Value expected to be one of the following: ' repmat('%s, ',1,length(expectedBorderCollapse)-1) '%s. The value was ''%s'' instead.'],expectedBorderCollapse{:},x)); 
checkDimensions = @(x) assert(isempty(x) || any(strcmpi(x,expectedDimensions)),sprintf(['Value expected to be one of the following: ' repmat('%s, ',1,length(expectedDimensions)-1) '%s. The value was ''%s'' instead.'],expectedDimensions{:},x)); 

addRequired(p,'Matrix',@(x) assert(isempty(x) || isnumeric(x) || iscell(x) || isstring(x),'Value must be a numeric, string, or cell array. Cells must contain numeric, character, or string values.'));
addParameter(p,'PageBackgroundColor',default.pagebackgroundcolor,checkColor1D)
addParameter(p,'Title',default.title,checkString1D)
addParameter(p,'TitleFontName',default.titlefontname,checkFonts)
addParameter(p,'TitleFontSize',default.titlefontsize,checkScalar)
addParameter(p,'TitleFontWeight',default.titlefontweight,checkWeights)
addParameter(p,'TitleFontStyle',default.titlefontstyle,checkStyle)
addParameter(p,'TitleFontColor',default.titlefontcolor,checkColor1D)
addParameter(p,'TitleTextAlign',default.titletextalign,checkAlignments)
addParameter(p,'TitleBackgroundColor',default.titlebackgroundcolor,checkColor1D)
addParameter(p,'LabelLayer',default.labellayer,checkLogicalScalar)
addParameter(p,'ColumnNames',default.columnnames,@(x) assert(isempty(x) || ((iscell(x) || isnumeric(x) || isstring(x)) && length(x)==size(matrix,2)),...
    'Value must be a [1xm] or [mx1] numeric, string, or cell array, where ''m'' is the number of columns in the data matrix. A cell array may contain numeric, character, or string values.'))
addParameter(p,'ColumnFontName',default.columnfontname,checkFonts)
addParameter(p,'ColumnFontSize',default.columnfontsize,checkScalar)
addParameter(p,'ColumnPadding',default.columnpadding,checkScalar)
addParameter(p,'ColumnFormat',default.columnformat,@(x) assert(isempty(x) || (ischar(x) && length(x)==size(matrix,2) && length(x) == numel(x)) || (iscell(x) && (length(x)==1 || (length(x)==size(matrix,2) && length(x) == numel(x)))),...
    'Value must be a valid sprintf output format, e.g. ''%.3f''. Valid input formats: 1D char array (applied to all column headings), [1x1] cell array (applied to all column headings), [1xm] or [mx1] cell array, where ''m'' is the number of columns in the data matrix (individual column heading formats).'))
addParameter(p,'ColumnFontWeight',default.columnfontweight,checkWeights)
addParameter(p,'ColumnFontStyle',default.columnfontstyle,checkStyle)
addParameter(p,'ColumnFontColor',default.columnfontcolor,@(x) assert(isempty(x) || (isnumeric(x) && (size(x,1)==1 || size(x,1)==size(matrix,2)) && size(x,2)==3 && all(x<=255 & x>=0,'all')),...
    'Value must be a [1x3] or [mx3] numeric array of RGB values, where ''m'' is the number of columns in the data matrix. Values must be scaled between 0 and 255.'))
addParameter(p,'ColumnMinimumWidth',default.columnminimumwidth,checkScalar)
addParameter(p,'ColumnTextAlign',default.columntextalign,checkAlignments)
addParameter(p,'ColumnBackgroundColor',default.columnbackgroundcolor,@(x) assert(isempty(x) || (isnumeric(x) && (size(x,1)==1 || size(x,1)==size(matrix,2)) && size(x,2)==3  && all(x<=255 & x>=0,'all')),...
    'Value must be a [1x3] or [mx3] numeric array of RGB values, where ''m'' is the number of columns in the data matrix. Values must be scaled between 0 and 255.'))
addParameter(p,'RowNames',default.rownames,@(x) assert(isempty(x) || ((iscell(x) || isnumeric(x) || isstring(x)) && length(x)==size(matrix,1)),...
    'Value must be a [1xn] or [nx1] numeric, string, or cell array, where ''n'' is the number of rows in the data matrix. A cell array may contain numeric, character, or string values.'))
addParameter(p,'RowFontName',default.rowfontname,checkFonts)
addParameter(p,'RowFontSize',default.rowfontsize,checkScalar)
addParameter(p,'RowPadding',default.rowpadding,checkScalar)
addParameter(p,'RowFormat',default.rowformat,@(x) assert(isempty(x) || (ischar(x) && length(x)==size(matrix,1) && length(x) == numel(x)) || (iscell(x) && (length(x)==1 || (length(x)==size(matrix,1) && length(x) == numel(x)))),...
    'Value must be a valid sprintf output format, e.g. ''%.3f''. Valid input formats: 1D char array (applied to all row headings), [1x1] cell array (applied to all row headings), [1xn] or [nx1] cell array, where ''n'' is the number of rows in the data matrix (individual row heading formats).'))
addParameter(p,'RowFontWeight',default.rowfontweight,checkWeights)
addParameter(p,'RowFontStyle',default.rowfontstyle,checkStyle)
addParameter(p,'RowFontColor',default.rowfontcolor,@(x) assert(isempty(x) || (isnumeric(x) && (size(x,1)==1 || size(x,1)==size(matrix,1)) && size(x,2)==3 && all(x<=255 & x>=0,'all')),...
    'Value must be a [1x3] or [nx3] numeric array of RGB values, where ''n'' is the number of rows in the data matrix. Values must be scaled between 0 and 255.'))
addParameter(p,'RowTextAlign',default.rowtextalign,checkAlignments)
addParameter(p,'RowBackgroundColor',default.rowbackgroundcolor,@(x) assert(isempty(x) || (isnumeric(x) && (size(x,1)==1 || size(x,1)==size(matrix,1)) && size(x,2)==3 && all(x<=255 & x>=0,'all')),...
    'Value must be a [1x3] or [nx3] numeric array or RGB values, where ''n'' is the number of rows in the data matrix. Values must be scaled between 0 and 255.'))
addParameter(p,'DataFontName',default.datafontname,checkFonts)
addParameter(p,'DataFontSize',default.datafontsize,checkScalar)
addParameter(p,'DataPadding',default.datapadding,checkScalar)
addParameter(p,'DataFormat',default.dataformat,@(x) assert(isempty(x) || (ischar(x) && length(x) == numel(x)) || (iscell(x) && (length(x)==1 || length(x)==size(matrix,2))),...
    'Value must be a valid sprintf output format, e.g. ''%.3f''. Valid input formats: 1D char array (applied to all data), [1x1] cell array (applied to all data), [1xm] or [mx1] cell array, where ''m'' is the number of columns in the data matrix (column-dependent data format).'))
addParameter(p,'DataFontWeight',default.datafontweight,checkWeights)
addParameter(p,'DataFontStyle',default.datafontstyle,checkStyle)
addParameter(p,'DataTextAlign',default.datatextalign,checkAlignments)
addParameter(p,'DataBackgroundColor',default.databackgroundcolor,@(x) assert(isempty(x) || ((isnumeric(matrix) || iscell(matrix)) && any(strcmpi(x,expectedColormaps))) || ((isnumeric(x) || iscell(matrix) || isstring(matrix)) && (all(size(x,1,2)==[1,3],'all') || all(size(x,1,2)==[size(matrix,1),3],'all') || all(size(x,1,2)==[size(matrix,2),3],'all') || all(size(x,1,2,3)==[size(matrix,1),size(matrix,2),3],'all') || all(size(x,1,2,3)==[size(matrix,1),size(matrix,2),3*size(matrix,3)],'all')) && all(x<=255 & x>=0,'all')),...
    'Value must be a [1x3], [nx3], [mx3], [nxmx3], or [n x m x 3k] numeric array of RGB values, where ''n'' is the number of rows in the data matrix, ''m'' is the number of columns in the data matrix, and ''k'' is the number of layers in the data matrix. Values must be scaled between 0 and 255. In the case of n = m, set the preferred color direction (''DataBackgroundColorDimension'') to ''row'' or ''column'' (default). Alternatively, if the matrix is numeric, or a cell array of numeric values, one can specify a supported  colormap (e.g. ''parula'') to scale the values to. The data will be scaled along the dimension specified by ''DataBackgroundColorDimension''. String arrays require an RGB color definition.'))
addParameter(p,'DataBackgroundColorDimension',default.databackgroundcolordimension,checkDimensions)
addParameter(p,'DataFontColor',default.datafontcolor,@(x) assert(isempty(x) || ((isnumeric(matrix) || iscell(matrix)) && any(strcmpi(x,expectedColormaps))) || ((isnumeric(x) || iscell(matrix) || isstring(matrix)) && (all(size(x,1,2)==[1,3],'all') || all(size(x,1,2)==[size(matrix,1),3],'all') || all(size(x,1,2)==[size(matrix,2),3],'all') || all(size(x,1,2,3)==[size(matrix,1),size(matrix,2),3],'all') || all(size(x,1,2,3)==[size(matrix,1),size(matrix,2),3*size(matrix,3)],'all')) && all(x<=255 & x>=0,'all')),...
    'Value must be a [1x3], [nx3], [mx3], [nxmx3], or [n x m x 3k] numeric array of RGB values, where ''n'' is the number of rows in the data matrix, ''m'' is the number of columns in the data matrix, and ''k'' is the number of layers in the data matrix. Values must be scaled between 0 and 255. In the case of n = m, set the preferred color direction (''DataBackgroundColorDimension'') to ''row'' or ''column'' (default). Alternatively, if the matrix is numeric, or a cell array of numeric values, one can specify a supported  colormap (e.g. ''parula'') to scale the values to. The data will be scaled along the dimension specified by ''DataBackgroundColorDimension''. String arrays require an RGB color definition.'))
addParameter(p,'DataFontColorDimension',default.datafontcolordimension,checkDimensions)
addParameter(p,'NegativeDataFontColor',default.negativedatafontcolor,checkColor1D)
addParameter(p,'DataEditable',default.dataeditable,checkLogical1D)
addParameter(p,'OmitValues',default.omitvalues,checkLogical1D)
addParameter(p,'Striping',default.striping,checkLogicalScalar)
addParameter(p,'StripingColors',default.stripingcolors,@(x) assert(isempty(x) || (isnumeric(x) && all(size(x)==[2,3]) && all(x<=255 & x>=0,'all')),...
    'Value must be a [2x3] numeric array of RGB values scaled between 0 and 255. The color defined in row one is applied to odd rows of the table. The color defined in row two is applied to even rows of the table.'))
addParameter(p,'BorderStyle',default.borderstyle,checkBorderStyles)
addParameter(p,'BorderWidth',default.borderwidth,checkScalar)
addParameter(p,'BorderColor',default.bordercolor,checkColor1D)
addParameter(p,'BorderCollapse',default.bordercollapse,checkBorderCollapse)
addParameter(p,'HighlightHover',default.highlighthover,checkLogicalScalar)
addParameter(p,'HoverColor',default.hovercolor,checkColor1D)
addParameter(p,'HighlightRows',default.highlightrows,@(x) assert(isempty(x) || (isnumeric(x) && length(x) == numel(x)),...
    'Value must be a 1D numeric array containing the row numbers to be highlighted. Negative numbers are ignored. Values are rounded to the nearest integer.'))
addParameter(p,'HighlightColor',default.highlightcolor,@(x) assert(isempty(x) || (isnumeric(x) && size(x,2)==3 && all(x<=255 & x>=0,'all')),...
    'Value must be a [1x3] or [nx3] numeric array of RGB values, where ''n'' is the number of values in the numeric array of ''HighlightRows''. Values must be scaled between 0 and 255.'))
addParameter(p,'ShowOutput',default.showoutput,checkLogicalScalar)
addParameter(p,'RecordRowClick',default.recordrowclick,checkLogicalScalar)
addParameter(p,'FileName',default.filename,checkString1D)
addParameter(p,'OutputToFile',default.outputtofile,checkLogicalScalar)
parse(p,matrix,varargin{:});

%======================================================
% Replace empty fields with their default values
%======================================================

param = p.Results;
fn = fieldnames(param);
for i = 1:length(fn)
    if isempty(param.(fn{i}))
        param.(fn{i}) = default.(lower(fn{i}));
    end
end

%======================================================
% Format parameters for consistency
% Check parameters and issue warnings/errors
%======================================================

matrixsize = size(param.Matrix);
param.Matrix = reshape(param.Matrix,matrixsize(1),matrixsize(2),prod(matrixsize(3:end)));

if param.Striping && isnumeric(param.Matrix) && ~isempty(param.NegativeDataFontColor)
    warning('HTMLtable: NegativeDataFontColor: Parameter is not used while row striping is active.');
end
if (iscell(param.Matrix) || isstring(param.Matrix)) && ~isempty(param.NegativeDataFontColor)
    warning('HTMLtable: NegativeDataFontColor: To apply, the parameter requires a numerical data array as the input.');
end
if ~isempty(param.HighlightRows) && (size(param.HighlightColor,1) ~= 1 && size(param.HighlightColor,1) ~= length(param.HighlightRows))
    error('HighlightColor: The number of rows of HighlightColor must be 1 or match the length of HighlightRows.')
end
param.HighlightRows = round(param.HighlightRows);
if param.Striping && ~isempty(param.DataBackgroundColor)
    warning('HTMLtable: DataBackgroundColor: Parameter is not used while row striping is active.');
end

ColormapMultiplier = 2;

if ~ischar(param.DataBackgroundColor)
    if size(param.Matrix,1) ~= size(param.Matrix,2)
        if isempty(param.DataBackgroundColor) || all(size(param.DataBackgroundColor,1,2)==[1,3],'all') || ...
                all(size(param.DataBackgroundColor,1,2,3)==[size(param.Matrix,1),size(param.Matrix,2),3],'all') || ...
                all(size(param.DataBackgroundColor,1,2,3)==[size(param.Matrix,1),size(param.Matrix,2),3*size(param.Matrix,3)],'all')
            param.DataBackgroundColorDimension = 'table';
        elseif all(size(param.DataBackgroundColor,1,2)==[size(param.Matrix,1),3],'all') 
            param.DataBackgroundColorDimension = 'row';
        elseif all(size(param.DataBackgroundColor,1,2)==[size(param.Matrix,2),3],'all') 
            param.DataBackgroundColorDimension = 'column';
        end
    else
        if isempty(param.DataBackgroundColor) || all(size(param.DataBackgroundColor,1,2)==[1,3],'all') || ...
                all(size(param.DataBackgroundColor,1,2,3)==[size(param.Matrix,1),size(param.Matrix,2),3],'all') || ...
                all(size(param.DataBackgroundColor,1,2,3)==[size(param.Matrix,1),size(param.Matrix,2),3*size(param.Matrix,3)],'all')
            param.DataBackgroundColorDimension = 'table';
        elseif all(size(param.DataBackgroundColor,1,2)==[size(param.Matrix,1),3],'all')
            if ~any(strcmpi(param.DataBackgroundColorDimension,{'row','column'}))
                param.DataBackgroundColorDimension = 'column';
            end
        end
    end
elseif ~param.Striping
    % Scale the values to the colormap along the specified dimension
    
    DBCfunc = str2func(param.DataBackgroundColor);
    
    if any(strcmpi(param.DataBackgroundColor,{'cividis','inferno','magma','plasma','twilight','viridis'}))
        try
            value = DBCfunc(1);
        catch
            error('HTMLtable: DataBackgroundColor: %s is not a built-in function. It can be downloaded on file exchange from several sources (e.g. Perceptually uniform colormaps by Ander Biguri).',param.DataBackgroundColor);
        end
    end
    
    switch param.DataBackgroundColorDimension
        case {'row'}
            nColorValues = ColormapMultiplier*size(param.Matrix,2);
            colors = DBCfunc(nColorValues);
            
            for ilayer = 1:size(param.Matrix,3)
                for i = 1:size(param.Matrix,1)
                    scaledValues = (param.Matrix(i,:,ilayer)-min(param.Matrix(i,:,ilayer),[],'all'))./(max(param.Matrix(i,:,ilayer),[],'all')-min(param.Matrix(i,:,ilayer),[],'all')) * nColorValues;
                    R = interp1(1:nColorValues,colors(:,1),scaledValues,'linear',colors(1,1));
                    G = interp1(1:nColorValues,colors(:,2),scaledValues,'linear',colors(2,1));
                    B = interp1(1:nColorValues,colors(:,3),scaledValues,'linear',colors(3,1));
                    param.DataBackgroundColor(i,1:size(param.Matrix,2),(1:3)+3*(ilayer-1)) = reshape([R,G,B]*255,1,size(param.Matrix,2),3);
                end
            end
            
        case {'column'}
            nColorValues = ColormapMultiplier*size(param.Matrix,1);
            colors = DBCfunc(nColorValues);
            
            for ilayer = 1:size(param.Matrix,3)
                for i = 1:size(param.Matrix,2)
                    scaledValues = (param.Matrix(:,i,ilayer)-min(param.Matrix(:,i,ilayer),[],'all'))./(max(param.Matrix(:,i,ilayer),[],'all')-min(param.Matrix(:,i,ilayer),[],'all')) * nColorValues;
                    R = interp1(1:nColorValues,colors(:,1),scaledValues,'linear',colors(1,1));
                    G = interp1(1:nColorValues,colors(:,2),scaledValues,'linear',colors(2,1));
                    B = interp1(1:nColorValues,colors(:,3),scaledValues,'linear',colors(3,1));
                    param.DataBackgroundColor(1:size(param.Matrix,1),i,(1:3)+3*(ilayer-1)) = reshape([R,G,B]*255,size(param.Matrix,1),1,3);
                end
            end
            
        case {'layer'}
            nColorValues = ColormapMultiplier*size(param.Matrix,1)*size(param.Matrix,2);
            colors = DBCfunc(nColorValues);
            param.DataBackgroundColor = [];
            
            for ilayer = 1:size(param.Matrix,3)
                scaledValues = (param.Matrix(:,:,ilayer)-min(param.Matrix(:,:,ilayer),[],'all'))./(max(param.Matrix(:,:,ilayer),[],'all')-min(param.Matrix(:,:,ilayer),[],'all')) * nColorValues;
                param.DataBackgroundColor(:,:,1+3*(ilayer-1)) = interp1(1:nColorValues,colors(:,1),scaledValues,'linear',colors(1,1))*255;
                param.DataBackgroundColor(:,:,2+3*(ilayer-1)) = interp1(1:nColorValues,colors(:,2),scaledValues,'linear',colors(2,1))*255;
                param.DataBackgroundColor(:,:,3+3*(ilayer-1)) = interp1(1:nColorValues,colors(:,3),scaledValues,'linear',colors(3,1))*255;
            end
            
        case {'table'}
            nColorValues = ColormapMultiplier*size(param.Matrix,1)*size(param.Matrix,2);
            colors = DBCfunc(nColorValues);
            param.DataBackgroundColor = [];
            scaledValues = (param.Matrix-min(param.Matrix,[],'all'))./(max(param.Matrix,[],'all')-min(param.Matrix,[],'all')) * nColorValues;
            
            for ilayer = 1:size(param.Matrix,3)
                param.DataBackgroundColor(:,:,1+3*(ilayer-1)) = interp1(1:nColorValues,colors(:,1),scaledValues(:,:,ilayer),'linear',colors(1,1))*255;
                param.DataBackgroundColor(:,:,2+3*(ilayer-1)) = interp1(1:nColorValues,colors(:,2),scaledValues(:,:,ilayer),'linear',colors(2,1))*255;
                param.DataBackgroundColor(:,:,3+3*(ilayer-1)) = interp1(1:nColorValues,colors(:,3),scaledValues(:,:,ilayer),'linear',colors(3,1))*255;
            end
    end
    param.DataBackgroundColorDimension = 'table';
end

if ~ischar(param.DataFontColor)
    if size(param.Matrix,1) ~= size(param.Matrix,2)
        if isempty(param.DataFontColor) || all(size(param.DataFontColor,1,2)==[1,3],'all') || all(size(param.DataFontColor,1,2,3)==[size(param.Matrix,1),size(param.Matrix,2),3],'all') || all(size(param.DataFontColor,1,2,3)==[size(param.Matrix,1),size(param.Matrix,2),3*size(param.Matrix,3)],'all')
            param.DataFontColorDimension = 'table';
        elseif all(size(param.DataFontColor,1,2)==[size(param.Matrix,1),3],'all') 
            param.DataFontColorDimension = 'row';
        elseif all(size(param.DataFontColor,1,2)==[size(param.Matrix,2),3],'all') 
            param.DataFontColorDimension = 'column';
        end
    else
        if isempty(param.DataFontColor) || all(size(param.DataFontColor,1,2)==[1,3],'all') || all(size(param.DataFontColor,1,2,3)==[size(param.Matrix,1),size(param.Matrix,2),3],'all') || all(size(param.DataFontColor,1,2,3)==[size(param.Matrix,1),size(param.Matrix,2),3*size(param.Matrix,3)],'all')
            param.DataFontColorDimension = 'table';
        elseif all(size(param.DataFontColor,1,2)==[size(param.Matrix,1),3],'all')
            if ~any(strcmpi(param.DataFontColorDimension,{'row','column'}))
                param.DataFontColorDimension = 'column';
            end
        end
    end
else
    % Scale the values to the colormap along the specified dimension
    
    DFCfunc = str2func(param.DataFontColor);
    
    if any(strcmpi(param.DataFontColor,{'cividis','inferno','magma','plasma','twilight','viridis'}))
        try
            value = DFCfunc(1);
        catch
            error('HTMLtable: DataFontColor: %s is not a built-in function. It can be downloaded on file exchange from several sources (e.g. Perceptually uniform colormaps by Ander Biguri).',param.DataFontColor);
        end
    end
    
    switch param.DataFontColorDimension
        case {'row'}
            nColorValues = ColormapMultiplier*size(param.Matrix,2);
            colors = DFCfunc(nColorValues);
            
            for ilayer = 1:size(param.Matrix,3)
                for i = 1:size(param.Matrix,1)
                    scaledValues = (param.Matrix(i,:,ilayer)-min(param.Matrix(i,:,ilayer),[],'all'))./(max(param.Matrix(i,:,ilayer),[],'all')-min(param.Matrix(i,:,ilayer),[],'all')) * nColorValues;
                    R = interp1(1:nColorValues,colors(:,1),scaledValues,'linear',colors(1,1));
                    G = interp1(1:nColorValues,colors(:,2),scaledValues,'linear',colors(2,1));
                    B = interp1(1:nColorValues,colors(:,3),scaledValues,'linear',colors(3,1));
                    param.DataFontColor(i,1:size(param.Matrix,2),(1:3)+3*(ilayer-1)) = reshape([R,G,B]*255,1,size(param.Matrix,2),3);
                end
            end
            
        case {'column'}
            nColorValues = ColormapMultiplier*size(param.Matrix,1);
            colors = DFCfunc(nColorValues);
            
            for ilayer = 1:size(param.Matrix,3)
                for i = 1:size(param.Matrix,2)
                    scaledValues = (param.Matrix(:,i,ilayer)-min(param.Matrix(:,i,ilayer),[],'all'))./(max(param.Matrix(:,i,ilayer),[],'all')-min(param.Matrix(:,i,ilayer),[],'all')) * nColorValues;
                    R = interp1(1:nColorValues,colors(:,1),scaledValues,'linear',colors(1,1));
                    G = interp1(1:nColorValues,colors(:,2),scaledValues,'linear',colors(2,1));
                    B = interp1(1:nColorValues,colors(:,3),scaledValues,'linear',colors(3,1));
                    param.DataFontColor(1:size(param.Matrix,1),i,(1:3)+3*(ilayer-1)) = reshape([R,G,B]*255,size(param.Matrix,1),1,3);
                end
            end
            
        case {'layer'}
            nColorValues = ColormapMultiplier*size(param.Matrix,1)*size(param.Matrix,2);
            colors = DFCfunc(nColorValues);
            param.DataFontColor = [];
            
            for ilayer = 1:size(param.Matrix,3)
                scaledValues = (param.Matrix(:,:,ilayer)-min(param.Matrix(:,:,ilayer),[],'all'))./(max(param.Matrix(:,:,ilayer),[],'all')-min(param.Matrix(:,:,ilayer),[],'all')) * nColorValues;
                param.DataFontColor(:,:,1+3*(ilayer-1)) = interp1(1:nColorValues,colors(:,1),scaledValues,'linear',colors(1,1))*255;
                param.DataFontColor(:,:,2+3*(ilayer-1)) = interp1(1:nColorValues,colors(:,2),scaledValues,'linear',colors(2,1))*255;
                param.DataFontColor(:,:,3+3*(ilayer-1)) = interp1(1:nColorValues,colors(:,3),scaledValues,'linear',colors(3,1))*255;
            end
            
        case {'table'}
            nColorValues = 10*size(param.Matrix,1)*size(param.Matrix,2);
            colors = DFCfunc(nColorValues);
            param.DataFontColor = [];
            scaledValues = (param.Matrix-min(param.Matrix,[],'all'))./(max(param.Matrix,[],'all')-min(param.Matrix,[],'all')) * nColorValues;
            
            for ilayer = 1:size(param.Matrix,3)
                param.DataFontColor(:,:,1+3*(ilayer-1)) = interp1(1:nColorValues,colors(:,1),scaledValues(:,:,ilayer),'linear',colors(1,1))*255;
                param.DataFontColor(:,:,2+3*(ilayer-1)) = interp1(1:nColorValues,colors(:,2),scaledValues(:,:,ilayer),'linear',colors(2,1))*255;
                param.DataFontColor(:,:,3+3*(ilayer-1)) = interp1(1:nColorValues,colors(:,3),scaledValues(:,:,ilayer),'linear',colors(3,1))*255;
            end
    end
    param.DataFontColorDimension = 'table';
end

% number of columns to be displayed (including header, if applicable)
ncols = size(param.Matrix,2) + double(~isempty(param.RowNames));
nHeaderRows = double(~isempty(param.Title)) + double(~isempty(param.ColumnNames)) + double(size(param.Matrix,3)>1);
matname = inputname(1);

if param.RecordRowClick
    rowclick = ' onclick="writeIndex(this.parentNode.rowIndex)"';
    onblur = ' onblur="writeValue(this.innerHTML)"';
else
    rowclick = '';
    onblur = '';
end

%======================================================
% Generate HTML
%======================================================
HTML = '<!DOCTYPE html>';

%--------------
% CSS
%--------------
HTML = [HTML '<STYLE>'];

HTML = [HTML,...
    sprintf('body{background-color: rgb(%d,%d,%d);}\n',param.PageBackgroundColor),...
    sprintf('Table,TR,TD,TH{border-style: %s; border-width: %dpx; border-color: rgb(%d,%d,%d); border-collapse: %s;}\n',...
        param.BorderStyle,param.BorderWidth,param.BorderColor(1,:),param.BorderCollapse),...
    sprintf('TD{font-family: %s; font-size: %dpt; font-weight: %s; font-style: %s; padding: %dpx; text-align: %s; min-width: %dpx;}\n',...
        param.DataFontName,param.DataFontSize,param.DataFontWeight,param.DataFontStyle,param.DataPadding,param.DataTextAlign,param.ColumnMinimumWidth),...
    sprintf('.title{color: rgb(%d,%d,%d); font-family: %s; font-size: %dpt; font-weight: %s; font-style: %s; text-align: %s; background-color: rgb(%d,%d,%d);}\n',...
        param.TitleFontColor(1,:),param.TitleFontName,param.TitleFontSize,param.TitleFontWeight,param.TitleFontStyle,param.TitleTextAlign,param.TitleBackgroundColor(1,:)),...
    sprintf('.layer{color: rgb(%d,%d,%d); background-color: rgb(%d,%d,%d);}\n',param.TitleFontColor(1,:),param.TitleBackgroundColor(1,:)),...
    sprintf('.emptycolumn{background-color: rgb(%d,%d,%d); min-width: %dpx;}\n',param.ColumnBackgroundColor(1,:),param.ColumnMinimumWidth)...
        ];
    
shortcutColumnClass = false;
if ~isempty(param.ColumnNames)
    for jj = 1:length(param.ColumnNames) 
        if size(param.ColumnFontColor,1) == 1
            CFC = param.ColumnFontColor(1,:);
        else
            CFC = param.ColumnFontColor(jj,:);
        end

        if size(param.ColumnBackgroundColor,1) == 1
            CBC = param.ColumnBackgroundColor(1,:);
        else
            CBC = param.ColumnBackgroundColor(jj,:);
        end

        if size(param.ColumnFontColor,1) == 1 && size(param.ColumnBackgroundColor,1) == 1
            HTML = [HTML sprintf('.columnLabel{color: rgb(%d,%d,%d); font-family: %s; font-size: %dpt; font-weight: %s; font-style: %s; padding: %dpx; text-align: %s; background-color: rgb(%d,%d,%d); min-width: %dpx;}\n',...
                    CFC,param.ColumnFontName,param.ColumnFontSize,param.ColumnFontWeight,param.ColumnFontStyle,param.ColumnPadding,param.ColumnTextAlign,CBC,param.ColumnMinimumWidth)];
            shortcutColumnClass = true;
            break
        else
            HTML = [HTML sprintf('.columnLabel%d{color: rgb(%d,%d,%d); font-family: %s; font-size: %dpt; font-weight: %s; font-style: %s; padding: %dpx; text-align: %s; background-color: rgb(%d,%d,%d); min-width: %dpx;}\n',...
                    jj,CFC,param.ColumnFontName,param.ColumnFontSize,param.ColumnFontWeight,param.ColumnFontStyle,param.ColumnPadding,param.ColumnTextAlign,CBC,param.ColumnMinimumWidth)];
        end
    end
end

shortcutRowClass = false;
if ~isempty(param.RowNames)
    for ii = 1:length(param.RowNames) 
        if size(param.RowFontColor,1) == 1
            RFC = param.RowFontColor(1,:);
        else
            RFC = param.RowFontColor(ii,:);
        end

        if size(param.RowBackgroundColor,1) == 1
            RBC = param.RowBackgroundColor(1,:);
        else
            RBC = param.RowBackgroundColor(ii,:);
        end

        if size(param.RowFontColor,1) == 1 && size(param.RowBackgroundColor,1) == 1
            HTML = [HTML sprintf('.rowLabel{color: rgb(%d,%d,%d); font-family: %s; font-size: %dpt; font-weight: %s; font-style: %s; padding: %dpx; text-align: %s; background-color: rgb(%d,%d,%d);}\n',...
                    RFC,param.RowFontName,param.RowFontSize,param.RowFontWeight,param.RowFontStyle,param.RowPadding,param.RowTextAlign,RBC)];
            shortcutRowClass = true;
            break
        else
            HTML = [HTML sprintf('.rowLabel%d{color: rgb(%d,%d,%d); font-family: %s; font-size: %dpt; font-weight: %s; font-style: %s; padding: %dpx; text-align: %s; background-color: rgb(%d,%d,%d);}\n',...
                    ii,RFC,param.RowFontName,param.RowFontSize,param.RowFontWeight,param.RowFontStyle,param.RowPadding,param.RowTextAlign,RBC)];
        end
    end
end

shortcutDataClass = false;
for ilayer = 1:size(param.Matrix,3)
    for jj = 1:size(param.Matrix,2)
        for ii = 1:size(param.Matrix,1)

            %DATA FONT COLOR
            if ~isempty(param.NegativeDataFontColor) && isnumeric(param.Matrix) && param.Matrix(ii,jj,ilayer) < 0 
                DFC = param.NegativeDataFontColor;
            else
                switch param.DataFontColorDimension
                    case {'row'}
                        DFC = param.DataFontColor(ii,:);
                    case {'column'}
                        DFC = param.DataFontColor(jj,:);
                    case {'table'}
                        if numel(param.DataFontColor) == 3
                            DFC = param.DataFontColor;
                        elseif size(param.DataFontColor,3) > 3
                            DFC = param.DataFontColor(ii,jj,(1:3)+3*(ilayer-1));
                        elseif size(param.DataFontColor,3) == 3
                            DFC = param.DataFontColor(ii,jj,:);
                        end
                end
            end

            if param.Striping
                if iscell(param.Matrix) || (isnumeric(param.Matrix) && numel(param.DataFontColor) == 3 && isempty(param.NegativeDataFontColor))
                    HTML = [HTML sprintf('.data{color: rgb(%d,%d,%d);}\n',DFC)];
                    shortcutDataClass = true;
                    break
                else
                    HTML = [HTML sprintf('.dataP%dC%dR%d{color: rgb(%d,%d,%d);}\n',ilayer,jj,ii,DFC)];
                end
            else
                
                %DATA BACKGROUND COLOR
                switch param.DataBackgroundColorDimension
                    case {'row'}
                        DBC = param.DataBackgroundColor(ii,:);
                    case {'column'}
                        DBC = param.DataBackgroundColor(jj,:);
                    case {'table'}
                        if numel(param.DataBackgroundColor) == 3
                            DBC = param.DataBackgroundColor;
                        elseif size(param.DataBackgroundColor,3) > 3
                            DBC = reshape(param.DataBackgroundColor(ii,jj,(1:3)+3*(ilayer-1)),[1,3,1]);
                        elseif size(param.DataBackgroundColor,3) == 3
                            DBC = reshape(param.DataBackgroundColor(ii,jj,:),[1,3,1]);
                        end
                end
                
                if numel(param.DataBackgroundColor) == 3 && numel(param.DataFontColor) == 3 && isempty(param.NegativeDataFontColor)
                    HTML = [HTML sprintf('.data{color: rgb(%d,%d,%d); background-color: rgb(%d,%d,%d);}\n',DFC,DBC)];
                    shortcutDataClass = true;
                    break
                else
                    HTML = [HTML sprintf('.dataP%dC%dR%d{color: rgb(%d,%d,%d); background-color: rgb(%d,%d,%d);}\n',ilayer,jj,ii,DFC,DBC)];
                end
            end
        end
        
        if shortcutDataClass
            break
        end
    end
    
    if shortcutDataClass
        break
    end
end

if ~isempty(param.HighlightRows)
    for i = 1:length(param.HighlightRows)
        if size(param.HighlightColor,1) == 1
            HTML = [HTML sprintf('.dataHighlight{background-color: rgb(%d,%d,%d);}\n',param.HighlightColor)];
        else
            HTML = [HTML sprintf('.dataHighlight%d{background-color: rgb(%d,%d,%d);}\n',param.HighlightRows(i),param.HighlightColor(i,:))];
        end
    end
end

if param.Striping
    HTML = [HTML sprintf('TR:nth-child(2n+%d){background-color: rgb(%d,%d,%d);}\n',nHeaderRows-1,param.StripingColors(1,:))];
	HTML = [HTML sprintf('TR:nth-child(2n+%d){background-color: rgb(%d,%d,%d);}\n',nHeaderRows,param.StripingColors(2,:))];
end

if param.HighlightHover
    HTML = [HTML sprintf('.highlight:hover{background-color: rgb(%d,%d,%d)!important;}\n',param.HoverColor(1,:))];
    HTML = [HTML sprintf('.highlight.dataHighlight:hover{background-color: rgb(%d,%d,%d)!important;}\n',param.HoverColor(1,:))];
end

HTML = [HTML '</STYLE>'];

% HTML table
HTML = [HTML '<TABLE id="table">'];

%--------------
% TITLE
%--------------
if ~isempty(param.Title)
    HTML = [HTML sprintf('<TR><TH class="title" COLSPAN=%g >%s</TH></TR>',ncols,param.Title)];
end

% cycle through layers (k) of the matrix
for ilayer = 1:size(param.Matrix,3)
    
    layer = param.Matrix(:,:,ilayer);
    
    if size(param.Matrix,3) > 1 && param.LabelLayer
        HTML = [HTML sprintf('<TR><TH class="layer" COLSPAN=%g>%s( : , : , %g )</TH></TR>',ncols,matname,ilayer)];
    end

    %--------------
    % COLUMN LABEL
    %--------------
    if ~isempty(param.ColumnNames)
        HTML = [HTML sprintf('<TR>')];
        if ~isempty(param.RowNames) && length(param.ColumnNames)<=size(layer,2)
            HTML = [HTML '<TD class="emptycolumn"></TD>'];
        end
        
        for jj = 1:length(param.ColumnNames) 
            
            if shortcutColumnClass
                classname = 'columnLabel';
            else
                classname = sprintf('columnLabel%d',jj);
            end
            
            if isempty(param.ColumnFormat)
                % Automatically discern data type
                if iscell(param.ColumnNames)
                    if ischar(param.ColumnNames{jj}) || isstring(param.ColumnNames{jj})
                        HTML = [HTML sprintf('<TD class="%s">%s</TD>',classname,param.ColumnNames{jj})];
                    elseif isnumeric(param.ColumnNames{jj})
                        HTML = [HTML sprintf('<TD class="%s">%g</TD>',classname,param.ColumnNames{jj})];
                    else
                        warning('Column names must be numeric or character values. ColumnName #%d has an unrecognized data type.',jj)
                    end
                elseif isstring(param.ColumnNames)
                    HTML = [HTML sprintf('<TD class="%s">%s</TD>',classname,param.ColumnNames(jj))];
                elseif isnumeric(param.ColumnNames)        
                    HTML = [HTML sprintf('<TD class="%s">%d</TD>',classname,param.ColumnNames(jj))];
                else
                    warning('Column names must be numeric or cell array. ColumnNames is an unsuppported array type.')
                end
            else            
                % Apply the specified column format regardless of data type
                if iscell(param.ColumnNames)
                    name = param.ColumnNames{jj};
                elseif isnumeric(param.ColumnNames) || isstring(param.ColumnNames)
                    name = param.ColumnNames(jj);
                else
                    warning('Column names must be numeric or character values. ColumnName #%d has an unrecognized data type.',jj)
                    name = ''; 
                end

                if ischar(param.ColumnFormat)
                    HTML = [HTML sprintf(['<TD class="%s">' param.ColumnFormat '</TD>'],classname,name)];
                elseif iscell(param.ColumnFormat) && length(param.ColumnFormat)==1
                    if isempty(param.ColumnFormat{1})
                        warning('Column formats must be character values. ColumnFormat #%d is empty.',jj)
                    end
                    HTML = [HTML sprintf(['<TD class="%s">' param.ColumnFormat{1} '</TD>'],classname,name)];
                elseif iscell(param.ColumnFormat) && length(param.ColumnFormat)==size(param.Matrix,2)
                    if isempty(param.ColumnFormat{jj})
                        warning('Column formats must be character values. ColumnFormat #%d is empty.',jj)
                    end
                    HTML = [HTML sprintf(['<TD class="%s">' param.ColumnFormat{jj} '</TD>'],classname,name)];
                else
                    warning('Column formats must be a character or cell array. ColumnFormat is an unsuppported array type.')
                end
            end
        end
        HTML = [HTML sprintf('</TR>\n')];
    end
    
    %--------------
    % ROW DATA
    %--------------
    
    for ii = 1:size(layer,1)
        
        HTML = [HTML '<TR class="highlight">'];  % new row of data

        %--------------
        % ROW LABEL
        %--------------

        if ~isempty(param.RowNames)
            
            if shortcutRowClass
                classname = 'rowLabel';
            else
                classname = sprintf('rowLabel%d',ii);
            end
            
            if isempty(param.RowFormat)
                % Automatically discern data type
                if iscell(param.RowNames)
                    if ischar(param.RowNames{ii}) || isstring(param.RowNames{ii})
                        HTML = [HTML sprintf('<TD class="%s">%s</TD>',classname,param.RowNames{ii})];
                    elseif isnumeric(param.RowNames{ii})
                        HTML = [HTML sprintf('<TD class="%s">%g</TD>',classname,param.RowNames{ii})];
                    else
                        warning('Row names must be numeric or character values. RowName #%d has an unrecognized data type.',ii)
                    end
                elseif isstring(param.RowNames)
                    HTML = [HTML sprintf('<TD class="%s">%s</TD>',classname,param.RowNames(ii))];
                elseif isnumeric(param.RowNames)    
                    HTML = [HTML sprintf('<TD class="%s">%d</TD>',classname,param.RowNames(ii))];
                else
                    warning('Row names must be numeric or cell array. RowNames is an unsuppported array type.')
                end
            else            
                % Apply the specified Row format regardless of data type
                if iscell(param.RowNames)
                    name = param.RowNames{ii};
                elseif isnumeric(param.RowNames) || isstring(param.RowNames)
                    name = param.RowNames(ii);
                else
                    warning('Row names must be numeric or character values. RowName #%d has an unrecognized data type.',ii)
                    name = ''; 
                end

                if ischar(param.RowFormat)
                    HTML = [HTML sprintf(['<TD class="%s">' param.RowFormat '</TD>'],classname,name)];
                elseif iscell(param.RowFormat) && length(param.RowFormat)==1
                    if isempty(param.RowFormat{1})
                        warning('Row formats must be character values. RowFormat #%d is empty.',ii)
                    end
                    HTML = [HTML sprintf(['<TD class="%s">' param.RowFormat{1} '</TD>'],classname,name)];
                elseif iscell(param.RowFormat) && length(param.RowFormat)==size(param.Matrix,1)
                    if isempty(param.RowFormat{ii})
                        warning('Row formats must be character values. RowFormat #%d is empty.',ii)
                    end
                    HTML = [HTML sprintf(['<TD class="%s">' param.RowFormat{ii} '</TD>'],classname,name)];
                else
                    warning('Row formats must be a character or cell array. RowFormat is an unsuppported array type.')
                end
            end
        end

        %--------------
        % CELL DATA
        %--------------
        for jj = 1:length(layer(ii,:)) % columns of data
            
            %DATA VALUE
            if iscell(layer) % if data is of type cell array
                
                % Extract data if single element cell array
                if iscell(layer{ii,jj}) && length(layer{ii,jj}) == 1
                    layer{ii,jj} = layer{ii,jj}{1};
                end
                
                celldata = layer{ii,jj};
                
            elseif isnumeric(layer) || isstring(layer)
                
                celldata = layer(ii,jj);
            end
            
            if shortcutDataClass
                classname = 'data';
            else
                classname = sprintf('dataP%dC%dR%d',ilayer,jj,ii);
            end
            
            if ~isempty(param.HighlightRows) && ~isempty(param.HighlightColor) && any(ii==param.HighlightRows,'all')
                if size(param.HighlightColor,1) == 1
                    classname = [classname ' dataHighlight'];
                else
                    classname = [classname sprintf(' dataHighlight%d',ii)];
                end
            end
            
            if param.HighlightHover
                classname = [classname ' highlight'];
            end
                
            %DATA FORMAT
            if isempty(param.DataFormat)
                % Automatically discern data type

                if ischar(celldata) || isstring(celldata)
                    DF = '%s';
                elseif isnumeric(celldata)
                    DF = '%g';
                else
                    warning('Matrix values must be numeric or character. Cell (%d,%d) has an unsupported data type.',ii,jj)
                    HTML = [HTML '<TD></TD>']; % add blank cell
                    continue
                end
            else
                % Apply the specified DataFormat to each column regardless of data type
                if ischar(param.DataFormat)
                    DF = param.DataFormat;
                elseif iscell(param.DataFormat) && length(param.DataFormat)==1
                    if isempty(param.DataFormat{1})
                        warning('Data formats must be character values. DataFormat is empty.')
                    end
                    DF = param.DataFormat{1};
                elseif iscell(param.DataFormat) && length(param.DataFormat)==size(param.Matrix,2)
                    if isempty(param.DataFormat{jj})
                        warning('Data formats must be character values. DataFormat #%d is empty.',jj)
                    end
                    DF = param.DataFormat{jj};
                else
                    warning('Data formats must be a character or cell array. DataFormat is an unsuppported array type.')
                end
            end

            if (isscalar(param.OmitValues) && param.OmitValues) || (length(param.OmitValues) == size(layer,2) && param.OmitValues(jj))
                if param.RecordRowClick
                    HTML = [HTML sprintf('<TD class="%s" onclick="writeIndex(this.parentNode.rowIndex)"></TD>',classname)];  
                else
                    HTML = [HTML sprintf('<TD class="%s"></TD>',classname)];  
                end
            else
                if (isscalar(param.DataEditable) && param.DataEditable) || (length(param.DataEditable) == size(layer,2) && param.DataEditable(jj))
                    editable = 'true';
                else
                    editable = 'false';
                end
                
                tabindex = ii+(jj-1)*size(layer,1)+(ilayer-1)*(size(layer,1)*size(layer,2));

                HTML = [HTML '<TD class="' classname '"' rowclick '><div contenteditable="' editable '" tabindex="' sprintf('%d',tabindex) '"' onblur '>' sprintf(DF,celldata) '</div></TD>']; %#ok<PFCEL>
            end
        end

        HTML = [HTML sprintf('\n</TR>')]; % close data row
    end
end
HTML = [HTML '</TABLE>'];  % close table

%--------------
% SCRIPTS
%--------------
if param.RecordRowClick
    script = '<script>'; 
    script = [script 'var oldrow = 0; var newrow = 0; var value = 0;'];
    script = [script 'function setup(Obj) {document.getElementById("table").addEventListener("click", function(event) {Obj.Data = [newrow,oldrow,value];});}'];
    script = [script 'function writeValue(v) {value = parseFloat(v);}'];
    script = [script 'function writeIndex(r) {oldrow = newrow;' sprintf('newrow = parseFloat(r-%d);',nHeaderRows-1) '}'];
    HTML = [HTML script '</script>']; 
end

%--------------
% DISPLAY OUTPUT, IF REQUESTED
%--------------
if param.ShowOutput
    web(['text://<html>' HTML '</html>'],'-notoolbar');
end

%--------------
% SAVE OUTPUT, IF REQUESTED
%--------------
if param.OutputToFile
    param.FileName = char(param.FileName);
    assert(~isempty(param.FileName),'HTMLtable: FileName: File name must be populated with a 1D char array or [1x1] string when outputing HTML to file.')
    if length(param.FileName) < 5 || ~strcmpi(param.FileName(end-4:end),'.html')
        filename = [param.FileName '.html'];
    else
        filename = param.FileName;
    end
    fid = fopen(filename,'w');
    fprintf(fid,'%s',HTML);
    fclose(fid);
end

end