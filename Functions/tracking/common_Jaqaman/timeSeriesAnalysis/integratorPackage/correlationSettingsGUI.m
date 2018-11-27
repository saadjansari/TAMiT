function varargout = correlationSettingsGUI(varargin)
% CORRELATIONSETTINGSGUI M-file for correlationSettingsGUI.fig
%      CORRELATIONSETTINGSGUI, by itself, creates a new CORRELATIONSETTINGSGUI or raises the existing
%      singleton*.
%
%      H = CORRELATIONSETTINGSGUI returns the handle to a new CORRELATIONSETTINGSGUI or the handle to
%      the existing singleton*.
%
%      CORRELATIONSETTINGSGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORRELATIONSETTINGSGUI.M with the given input arguments.
%
%      CORRELATIONSETTINGSGUI('Property','Value',...) creates a new CORRELATIONSETTINGSGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before correlationSettingsGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to correlationSettingsGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help correlationSettingsGUI

% Last Modified by GUIDE v2.5 01-Feb-2012 17:51:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @correlationSettingsGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @correlationSettingsGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before correlationSettingsGUI is made visible.
function correlationSettingsGUI_OpeningFcn(hObject, eventdata, handles, varargin)

set(handles.text_copyright, 'String', getLCCBCopyright());


userData = get(handles.figure1, 'UserData');
if isempty(userData), userData = struct(); end
handles.output = hObject;

% Get main figure handle and process id
userData.iTool = varargin{4};
userData.mainFig = varargin{3}.figure1;
userData_main = get(userData.mainFig, 'UserData');

% Parameter Setup
processingTools = userData_main.tools(userData.iTool);
set(handles.edit_nBoot, 'String', processingTools.parameters.nBoot)
set(handles.edit_alpha, 'String', processingTools.parameters.alpha);

% Get icon infomation
userData.questIconData = userData_main.questIconData;
userData.colormap = userData_main.colormap;

% ----------------------Set up help icon------------------------

% Set up help icon
set(hObject,'colormap',userData.colormap);
% Set up package help. Package icon is tagged as '0'
set(handles.figure1,'CurrentAxes',handles.axes_help);Img = image(userData.questIconData); 
set(gca, 'XLim',get(Img,'XData'),'YLim',get(Img,'YData'),...
    'visible','off','YDir','reverse');
set(Img,'ButtonDownFcn',@icon_ButtonDownFcn,...
    'UserData', struct('class', mfilename));

set(handles.figure1, 'UserData', userData)
% Update handles structure
guidata(hObject, handles);



% UIWAIT makes correlationSettingsGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = correlationSettingsGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1)


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

userData=get(handles.figure1,'UserData');
if isempty(userData), userData = struct(); end

nBoot=str2double(get(handles.edit_nBoot,'String'));
if isnan(nBoot) || nBoot<=0
    errordlg(sprintf('Invalid value for the %s.',get(handles.text_nBoot,'String')),'Setting Error','modal')
    return;
end

alpha=str2double(get(handles.edit_alpha,'String'));
if isnan(alpha) || alpha<0 || alpha >1
    errordlg(sprintf('Invalid value for the %s.',get(handles.text_alpha,'String')),'Setting Error','modal')
    return;
end

userData_main = get(userData.mainFig,'UserData');
userData_main.tools(userData.iTool).parameters.nBoot=nBoot;
userData_main.tools(userData.iTool).parameters.alpha=alpha;
set(userData.mainFig,'UserData',userData_main);

delete(handles.figure1);


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
if strcmp(eventdata.Key, 'return')
    pushbutton_done_Callback(handles.pushbutton_done, [], handles);
end
