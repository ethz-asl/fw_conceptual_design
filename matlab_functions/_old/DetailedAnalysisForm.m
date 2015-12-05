function varargout = DetailedAnalysisForm(varargin)
% DETAILEDANALYSISFORM MATLAB code for DetailedAnalysisForm.fig
%      DETAILEDANALYSISFORM, by itself, creates a new DETAILEDANALYSISFORM or raises the existing
%      singleton*.
%
%      H = DETAILEDANALYSISFORM returns the handle to a new DETAILEDANALYSISFORM or the handle to
%      the existing singleton*.
%
%      DETAILEDANALYSISFORM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DETAILEDANALYSISFORM.M with the given input arguments.
%
%      DETAILEDANALYSISFORM('Property','Value',...) creates a new DETAILEDANALYSISFORM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DetailedAnalysisForm_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DetailedAnalysisForm_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DetailedAnalysisForm

% Last Modified by GUIDE v2.5 24-May-2012 16:09:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DetailedAnalysisForm_OpeningFcn, ...
                   'gui_OutputFcn',  @DetailedAnalysisForm_OutputFcn, ...
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

% --- Executes just before DetailedAnalysisForm is made visible.
function DetailedAnalysisForm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DetailedAnalysisForm (see VARARGIN)

% Choose default command line output for DetailedAnalysisForm
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using DetailedAnalysisForm.
%--------------------------------------------------------------------------
% INITIALIZE THE GUI HERE
%--------------------------------------------------------------------------
%cla;

%Get data from varargin
AR_array=varargin{1}{1};
b_array=varargin{1}{2};
m_bat_array=varargin{1}{3};
global flightdata2_array;
global performance2;
flightdata2_array=varargin{1}{4};
performance2.t_excess=varargin{1}{5};
performance2.t_endurance=varargin{1}{6};
performance2.t_chargemargin=varargin{1}{7};

%PopupMenus
tmp=0;
for i=1:length(AR_array); tmp(i)=AR_array(i); end
set(handles.PM_AR, 'String',tmp);
set(handles.PM_AR, 'Value',1);
tmp=0;
for i=1:length(b_array); tmp(i)=b_array(i); end
set(handles.PM_b, 'String',tmp);
set(handles.PM_b, 'Value',1);
tmp=0;
for i=1:length(m_bat_array); tmp(i)=m_bat_array(i); end
set(handles.PM_mbat, 'String',tmp);
set(handles.PM_mbat, 'Value',1);

%Update with first dataset
PB_Update_Callback(hObject, eventdata, handles);

% UIWAIT makes DetailedAnalysisForm wait for user response (see UIRESUME)
 %uiwait(handles.figure1);




% --- Outputs from this function are returned to the command line.
function varargout = DetailedAnalysisForm_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in PB_Update.
function PB_Update_Callback(hObject, eventdata, handles)
% hObject    handle to PB_Update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%figure(handles.figure1);

cla(handles.axes1);
cla(handles.axes2);
cla(handles.axes3);
cla(handles.axes4);

%Get selected configuration
ARidx = get(handles.PM_AR, 'Value');
bidx = get(handles.PM_b, 'Value');
mbatidx = get(handles.PM_mbat, 'Value');

global flightdata2_array;
global performance2;

%-------------------------------------------------
% OUTPUT PERFORMANCE DATA
%-------------------------------------------------
tmp1=['b: ' num2str(flightdata2_array(bidx,mbatidx,ARidx).b) 'm'];
tmp2=['mbat: ' num2str(flightdata2_array(bidx,mbatidx,ARidx).m_bat) 'kg'];
tmp3=['AR: ' num2str(flightdata2_array(bidx,mbatidx,ARidx).AR)];

tmp4=['turbulence: ' num2str(flightdata2_array(bidx,mbatidx,ARidx).turbulence)];
tmp5=['clearness: ' num2str(flightdata2_array(bidx,mbatidx,ARidx).clearness)];

tmp6=['t_excess: ' num2str(performance2.t_excess(bidx,mbatidx,ARidx))];
tmp7=['t_endurance: ' num2str(performance2.t_endurance(bidx,mbatidx,ARidx))];
tmp8=['t_chargemargin: ' num2str(performance2.t_chargemargin(bidx,mbatidx,ARidx))];

%check if we have data to draw!!
if(isempty(flightdata2_array(bidx,mbatidx,ARidx).t_array))
    set(handles.E_Performance,'String',{tmp1 tmp2 tmp3 ' ' tmp4 tmp5 ' ' tmp6 tmp7 tmp8});
    return;
end

deltat=flightdata2_array(bidx,mbatidx,ARidx).t_array(2)-flightdata2_array(bidx,mbatidx,ARidx).t_array(1);
overallt=flightdata2_array(bidx,mbatidx,ARidx).t_array(length(flightdata2_array(bidx,mbatidx,ARidx).t_array))-flightdata2_array(bidx,mbatidx,ARidx).t_array(1);
distancetravelled=0;
for i=1:length(flightdata2_array(bidx,mbatidx,ARidx).t_array)
    distancetravelled=distancetravelled+flightdata2_array(bidx,mbatidx,ARidx).v_array(i)*deltat;
end
tmp9=['distance: ' num2str(distancetravelled/1000,'%.2f') 'km/' num2str(overallt/3600.0,'%.1f') 'h'];

%Find night duration
idx_sunset=find(flightdata2_array(bidx,mbatidx,ARidx).irr_array<=0.1,1,'first');
idx_sunrise=idx_sunset+find(flightdata2_array(bidx,mbatidx,ARidx).irr_array(idx_sunset:end)>0,1,'first');
if idx_sunrise>numel(flightdata2_array(bidx,mbatidx,ARidx).t_array) 
    idx_sunrise=idx_sunset;
end
t_sunset=flightdata2_array(bidx,mbatidx,ARidx).t_array(idx_sunset);
t_sunrise=flightdata2_array(bidx,mbatidx,ARidx).t_array(idx_sunrise);
tmp10=['t_sunrise: ' num2str(mod(t_sunrise,86400)/3600) 'hrs'];
tmp11=['t_sunset: ' num2str(mod(t_sunset,86400)/3600) 'hrs'];
tmp12=['T_night: ' num2str((t_sunrise-t_sunset)/3600) 'hrs'];

set(handles.E_Performance,'String',{tmp1 tmp2 tmp3 ' ' tmp4 tmp5 ' ' tmp6 tmp7 tmp8 ' ' tmp9 ' ' tmp10 tmp11 tmp12});

%-------------------------------------------------
% DRAW FLIGHTDATA vs TIME
%-------------------------------------------------

%figure 1
% [handles.axyy NULL NULL]= plotyy(flightdata2_array(bidx,mbatidx,ARidx).t_array/3600,flightdata2_array(bidx,mbatidx,ARidx).h_array,...
%        flightdata2_array(bidx,mbatidx,ARidx).t_array/3600,flightdata2_array(bidx,mbatidx,ARidx).v_array);
% legend('h[m]','v[m/ss]');
plot(handles.axes1, flightdata2_array(bidx,mbatidx,ARidx).t_array/3600,flightdata2_array(bidx,mbatidx,ARidx).h_array);
legend(handles.axes1, 'h[m]');

%figure2
dh_array=zeros(size(flightdata2_array(bidx,mbatidx,ARidx).h_array));
dh_array(1)=0;
for i=2:size(flightdata2_array(bidx,mbatidx,ARidx).h_array,2) 
    dh_array(i)=flightdata2_array(bidx,mbatidx,ARidx).h_array(i)-flightdata2_array(bidx,mbatidx,ARidx).h_array(i-1); 
end

plot(handles.axes2, ...
     flightdata2_array(bidx,mbatidx,ARidx).t_array/3600,flightdata2_array(bidx,mbatidx,ARidx).Re_array/100,...
     flightdata2_array(bidx,mbatidx,ARidx).t_array/3600,flightdata2_array(bidx,mbatidx,ARidx).CL_array*1000,...
     flightdata2_array(bidx,mbatidx,ARidx).t_array/3600,flightdata2_array(bidx,mbatidx,ARidx).CD_array*10000,...
     flightdata2_array(bidx,mbatidx,ARidx).t_array/3600,dh_array*10);
legend(handles.axes2,'Re/100','Cl*1000','Cd*10000','dh*10');

%figure 3
plot(handles.axes3,...
    flightdata2_array(bidx,mbatidx,ARidx).t_array/3600,flightdata2_array(bidx,mbatidx,ARidx).bat_array/3600);
legend(handles.axes3,'E_{Bat}[Wh]');
ylim(handles.axes3,[0 max(flightdata2_array(bidx,mbatidx,ARidx).bat_array/3600)*1.1]);

%figure 4
plot(handles.axes4,...
    flightdata2_array(bidx,mbatidx,ARidx).t_array/3600,flightdata2_array(bidx,mbatidx,ARidx).P_solar_array,...
    flightdata2_array(bidx,mbatidx,ARidx).t_array/3600,flightdata2_array(bidx,mbatidx,ARidx).P_elec_tot_array,...
    flightdata2_array(bidx,mbatidx,ARidx).t_array/3600,flightdata2_array(bidx,mbatidx,ARidx).P_prop_array);
    %flightdata2_array(bidx,mbatidx,ARidx).t_array/3600,P_0);
    %ylabel(handles.axes4,'Irradiance [W/m^2]');
    xlabel(handles.axes4,'Time of Day [h]');
    legend(handles.axes4,'P_{Solar}[W]','P_{electot}[W]','P_{prop}[W]');%,'P_0[W]');
        
linkaxes([handles.axes1 handles.axes2 handles.axes3 handles.axes4],'x');
% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in PM_AR.
function PM_AR_Callback(hObject, eventdata, handles)
% hObject    handle to PM_AR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns PM_AR contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PM_AR


% --- Executes during object creation, after setting all properties.
function PM_AR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PM_AR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in PM_mbat.
function PM_mbat_Callback(hObject, eventdata, handles)
% hObject    handle to PM_mbat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PM_mbat contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PM_mbat


% --- Executes during object creation, after setting all properties.
function PM_mbat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PM_mbat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in PM_b.
function PM_b_Callback(hObject, eventdata, handles)
% hObject    handle to PM_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PM_b contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PM_b


% --- Executes during object creation, after setting all properties.
function PM_b_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PM_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when uipanel1 is resized.
function uipanel1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function E_Performance_Callback(hObject, eventdata, handles)
% hObject    handle to E_Performance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of E_Performance as text
%        str2double(get(hObject,'String')) returns contents of E_Performance as a double


% --- Executes during object creation, after setting all properties.
function E_Performance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E_Performance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PB_ExportPerfData.
function PB_ExportPerfData_Callback(hObject, eventdata, handles)
% hObject    handle to PB_ExportPerfData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global flightdata2_array; %Indices b,mbat,AR
global performance2;      %same

exportdata=zeros(numel(performance2.t_excess),6);
dim1=size(performance2.t_excess,1);
dim2=size(performance2.t_excess,2);
dim3=size(performance2.t_excess,3);

for i=1:dim1
    for j=1:dim2
        for k=1:dim3
        exportdata((i-1)*dim2*dim3+(j-1)*dim3+k,1)=flightdata2_array(i,j,k).b;
        exportdata((i-1)*dim2*dim3+(j-1)*dim3+k,2)=flightdata2_array(i,j,k).m_bat;
        exportdata((i-1)*dim2*dim3+(j-1)*dim3+k,3)=flightdata2_array(i,j,k).AR;
        exportdata((i-1)*dim2*dim3+(j-1)*dim3+k,4)=performance2.t_excess(i,j,k);
        exportdata((i-1)*dim2*dim3+(j-1)*dim3+k,5)=performance2.t_endurance(i,j,k);
        exportdata((i-1)*dim2*dim3+(j-1)*dim3+k,6)=performance2.t_chargemargin(i,j,k);
        end
    end
end
xlswrite('results\Output.xls',exportdata);
