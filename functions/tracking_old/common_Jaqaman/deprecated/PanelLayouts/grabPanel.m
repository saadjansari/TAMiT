function fig = grabPanel()
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
% This problem is solved by saving the output as a FIG-file.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.
% 
% NOTE: certain newer features in MATLAB may not have been saved in this
% M-file due to limitations of this format, which has been superseded by
% FIG-files.  Figures which have been annotated using the plot editor tools
% are incompatible with the M-file/MAT-file format, and should be saved as
% FIG-files.

load grabPanel
fprintf(2,['Warning: ''' mfilename ''' is deprecated and should no longer be used.\n']);

h0 = figure('Color',[0.8 0.8 0.8], ...
	'Colormap',[], ...
	'CreateFcn','uiAcqPanelCreateFnc', ...
	'DeleteFcn','uiAcqPanelDeleteFnc', ...
	'FileName','C:\usr\thomann\matlab\gui\PanelLayouts\grabPanel.m', ...
	'MenuBar','none', ...
	'MinColormap',0, ...
	'Name','Acquisition Panel', ...
	'NumberTitle','off', ...
	'PaperPosition',[18 180 576 432], ...
	'PaperUnits','points', ...
	'Position',[699 230 315 325], ...
	'Tag','GRABWINDOW', ...
	'ToolBar','none', ...
	'UserData',mat0);
h1 = uimenu('Parent',h0, ...
	'Callback','uiAcqMenuFileFnc', ...
	'Label','&File', ...
	'Tag','UIACQMENU_FILE');
h2 = uimenu('Parent',h1, ...
	'Callback','uiAcqMenuLoadFnc', ...
	'Label','&Load ...', ...
	'Tag','UIACQMENU_FILE_LOAD');
h2 = uimenu('Parent',h1, ...
	'Label','Load Stack', ...
	'Tag','UIACQMENU_FILE_LOADSTACK');
h3 = uimenu('Parent',h2, ...
	'Callback','uiAcqMenuLoadFullStackFnc', ...
	'Label','Full ...', ...
	'Tag','UIACQMENU_FILE_LOADSTACK_FULL');
h3 = uimenu('Parent',h2, ...
	'Callback','uiAcqMenuLoadClippedStackFnc', ...
	'Label','Clipped ...', ...
	'Tag','UIACQMENU_FILE_LOADSTACK_CLIPPED');
h2 = uimenu('Parent',h1, ...
	'Callback','uiAcqMenuLoadMpegFnc', ...
	'Label','Load MPEG ...', ...
	'Tag','UIACQMENU_FILE_LOADMPEG');
h2 = uimenu('Parent',h1, ...
	'Callback','uiAcqMenuSaveFnc', ...
	'Enable','off', ...
	'Label','&Save ...', ...
	'Tag','UIACQMENU_FILE_SAVE');
h2 = uimenu('Parent',h1, ...
	'Callback','uiAcqMenuSavempegFnc', ...
	'Enable','off', ...
	'Label','Save MPEG ...', ...
	'Tag','UIACQMENU_FILE_SAVEMPEG');
h1 = uimenu('Parent',h0, ...
	'Callback','uiAcqMenuFrameFnc', ...
	'Label','Fra&me', ...
	'Tag','UIACQMENU_FRAME');
h2 = uimenu('Parent',h1, ...
	'Callback','uiAcqMenuFullFnc', ...
	'Checked','on', ...
	'Label','Full', ...
	'Tag','UIACQMENU_FRAME_FULL');
h2 = uimenu('Parent',h1, ...
	'Callback','uiAcqMenuSelectFnc', ...
	'Label','Select', ...
	'Tag','UIACQMENU_FRAME_SELECT');
h3 = uimenu('Parent',h2, ...
	'Callback','uiAcqMenuFinteractiveFnc', ...
	'Enable','off', ...
	'Label','Interactive', ...
	'Tag','UIACQMENU_FRAME_SELECT_INTERACTIVE');
h3 = uimenu('Parent',h2, ...
	'Callback','uiAcqMenuFpanelFnc', ...
	'Label','Panel', ...
	'Tag','UIACQMENU_FRAME_SELECT_PANEL');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'CreateFcn','uiAcqDirEditCreateFnc', ...
	'Enable','off', ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[16.5 183.75 165 15], ...
	'String','c:\usr\thomann\matlab\gui\panellayouts', ...
	'Style','edit', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'CreateFcn','uiAcqFileEditCreateFnc', ...
	'Enable','off', ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[16.5 147.75 165 15], ...
	'Style','edit', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','uiAcqGrabSaveButtonPressFnc', ...
	'CreateFcn','uiAcqGrabSaveButtonCreateFnc', ...
	'Enable','off', ...
	'ListboxTop',0, ...
	'Position',[114 127.5 60 15], ...
	'String','Grab to File', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','uiAcqGrabButtonPressFnc', ...
	'CreateFcn','uiAcqGrabButtonCreateFnc', ...
	'Enable','off', ...
	'ListboxTop',0, ...
	'Position',[16.5 126.75 60 15], ...
	'String','Grab to View', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','uiAcqGrabStackButtonPressFnc', ...
	'CreateFcn','uiAcqGrabStackButtonCreateFnc', ...
	'Enable','off', ...
	'ListboxTop',0, ...
	'Position',[16.5 106.5 60 15], ...
	'String','Grab Stack', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','uiAcqMovieButtonPressFnc', ...
	'CreateFcn','uiAcqMovieButtonCreateFnc', ...
	'Enable','off', ...
	'ListboxTop',0, ...
	'Position',[16.5 84 60 15], ...
	'String','Create Movie', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'Callback','uiAcqFileButtonPressFnc', ...
	'CreateFcn','uiAcqFileButtonCreateFnc', ...
	'Enable','off', ...
	'ListboxTop',0, ...
	'Position',[196.5 147.75 30 15], ...
	'String','...', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'Callback','uiAcqStackSliderFnc', ...
	'Enable','off', ...
	'ListboxTop',0, ...
	'Max',50, ...
	'Min',30, ...
	'Position',[16.5 54 165 15], ...
	'SliderStep',[0.05 1], ...
	'Style','slider', ...
	'Tag','StaticText1', ...
	'Value',30);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'CreateFcn','uiAcqTimeEditCreateFnc', ...
	'Enable','off', ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[114 105 30 15], ...
	'String','10', ...
	'Style','edit', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback','uiAcqStackListFnc', ...
	'Enable','off', ...
	'Position',[196.5 52.5 30 30], ...
	'Style','listbox', ...
	'Tag','StaticText1', ...
	'Value',1);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'CreateFcn','uiAcqFpsEditCreateFnc', ...
	'Enable','off', ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[196.5 105 30 15], ...
	'String','5', ...
	'Style','edit', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'ListboxTop',0, ...
	'Position',[87.75 106.5 21 12], ...
	'String','[sec]:', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[151.5 106.5 33.75 12], ...
	'String','[frm/sec]:', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'ListboxTop',0, ...
	'Position',[185.25 126.75 11.25 15], ...
	'Style','checkbox', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[200.25 129 24.75 11.25], ...
	'String','Show', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'CreateFcn','uiAcqRepeatEditCreateFnc', ...
	'Enable','off', ...
	'HorizontalAlignment','right', ...
	'ListboxTop',0, ...
	'Position',[114 84 29.25 15], ...
	'String','1', ...
	'Style','edit', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'HorizontalAlignment','left', ...
	'ListboxTop',0, ...
	'Position',[151.5 86.25 30 11.25], ...
	'String','times', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Enable','off', ...
	'Position',[106.5 221.25 120 15], ...
	'String',' N/A', ...
	'Style','listbox', ...
	'Tag','StaticText1', ...
	'Value',1);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'ListboxTop',0, ...
	'Position',[12.75 222.75 82.5 11.25], ...
	'String','Select Framegrabber:', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'ListboxTop',0, ...
	'Position',[16.5 202.5 75 11.25], ...
	'String','Acquisition Directory:', ...
	'Style','text', ...
	'Tag','StaticText2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'ListboxTop',0, ...
	'Position',[12.75 165 63.75 11.25], ...
	'String','Acquisition File:', ...
	'Style','text', ...
	'Tag','StaticText3');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[129.75 19.5 37.5 18.75], ...
	'String','100', ...
	'Style','edit', ...
	'Tag','SHUTTER_TIME');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.8 0.8 0.8], ...
	'ListboxTop',0, ...
	'Position',[18.75 18.75 93.75 18.75], ...
	'String','Shutter Time in microsecs:', ...
	'Style','text', ...
	'Tag','StaticText4');
if nargout > 0, fig = h0; end