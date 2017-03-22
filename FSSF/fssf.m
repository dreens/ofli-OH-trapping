function varargout = fssf(varargin)
% FSSF Application M-file for fssf.fig
%    FIG = FSSF launch fssf GUI.
%    FSSF('callback_name', ...) invoke the named callback.
%    Graphical User Interface for the FSSF Kernel
%    (c) Miltiades Salvanos 2003
%    Last Modified by GUIDE v2.0 23-May-2003 20:00:02

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end
    
    vars = evalin('base','who'); % get the list of variables from the current workspace
    set(handles.dataselector,'String',vars) % populate the dataselector with this list
    
    set(handles.text2,'Enable','off') %deactivate further steps help
    set(handles.text3,'Enable','off')
    set(handles.text4,'Enable','off')
    set(handles.maxdev,'Enable','off'); % deactivate specific controls
    set(handles.maxdevpc,'Enable','off');
    set(handles.allowselected,'Enable','off');
    set(handles.allowany,'Enable','off');
    set(handles.plot,'Enable','off');
    set(handles.accept,'Enable','off');
    set(handles.exportdxf,'Enable','off');
    set(handles.bcx0f,'Value',1); % set default parameters
    set(handles.bcxmf,'Value',1);
    set(handles.bcy0f,'Value',1);
    set(handles.bcymf,'Value',1);
    set(handles.orderselector,'Value',1);
    set(handles.seriesselector,'Value',1);
    set(handles.allowselected,'Value',1);
    set(handles.allowany,'Value',0);
    handles.maxdevvalue=0;
    handles.order=3; handles.series='SIN*SIN';
    guidata(fig, handles);
 
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
		end
    catch
        disp(lasterr);
	end

end


%| Additional FUNCTIONS

function dataset = get_var_names(handles) % subfunction to obtain dataset selection  
list_entries = get(handles.dataselector,'String');
index_selected = get(handles.dataselector,'Value');
dataset = list_entries{index_selected};

function order = get_order_selection(handles) % subfunction to obtain order selection  
list_entries = get(handles.orderselector,'String');
index_selected = get(handles.orderselector,'Value');
order= list_entries{index_selected};

function series= get_series_selection(handles) % subfunction to obtain series selection  
list_entries = get(handles.seriesselector,'String');
index_selected = get(handles.seriesselector,'Value');
series = list_entries{index_selected};

function coefficients=present_coefficients(coeff)
size=(size(coeff));
coeff_size=(size(1))^0.5;
coefficients=[]; horiz=1; vert=1;
for vert=1:coeff_size;
    for horiz=1:coeff_size
       coefficients(horiz,vert)=coeff(coeff_size*(horiz-1)+vert);  %row index is for x-direction, column index is for y-direction
       end
   end

function plot_results(x,y,z,zfitted) % subfunction to mesh and plot results
figure;
xline=linspace(min(x),max(x),50);
yline=linspace(min(y),max(y),50*floor(abs((max(x)-min(x))/(max(y)-min(y)))));
[X,Y]=meshgrid(xline,yline);
Z=griddata(x,y,z,X,Y,'v4');
Zfitted=griddata(x,y,zfitted,X,Y,'v4');
DevColour=(Z-Zfitted);
surface(X,Y,Zfitted,DevColour,'EdgeColor','white','EdgeAlpha',0.5,'FaceColor','interp');
axis on; axis vis3d; axis equal; axis tight; grid on; colorbar('horiz'); view(3);
title('Fourier Series fitted Surface-Colour by Deviation');
xlabel('a  (x)'); ylabel('b  (y)'); zlabel('w  (z)')

function update_dev_graph(handles) % subfunction to update the maximum deviation graph
x=handles.order;  y=(100*handles.Max_Dev)./max(abs(handles.dataset(:,3)));
plot(x,y,'r',x,y,'ok');
set(gca,'Color',[1,0.99,0.97]);
axis([0 20 0 1.5*max(y)]); grid on;

%| Object CALLBACKS:

% --------------------------------------------------------------------
function varargout = close_Callback(h, eventdata, handles, varargin)
close(handles.window)



% --------------------------------------------------------------------
function varargout = accept_Callback(h, eventdata, handles, varargin)
csvwrite([handles.dataset_name '.dat'],handles.coefficients);
evalin('base',[handles.dataset_name,'_coefficients=',mat2str(handles.coefficients)]);


% --------------------------------------------------------------------
function varargout = plot_Callback(h, eventdata, handles, varargin)
plot_results(handles.dataset(:,1),handles.dataset(:,2),handles.dataset(:,3),handles.zfitted);


% --------------------------------------------------------------------
function varargout = fit_Callback(h, eventdata, handles, varargin)
handles.bcx_0=2-get(handles.bcx0s,'Value'); % 1 for supported, 2 for free
handles.bcx_end=2-get(handles.bcx0s,'Value');
handles.bcy_0=2-get(handles.bcx0s,'Value');
handles.bcy_end=2-get(handles.bcx0s,'Value');
if (get(handles.optimization,'Value'))==0
    output=feval(@kernel,handles.dataset,handles.order,handles.series,handles.bcx_0,handles.bcx_end,handles.bcy_0,handles.bcy_end);
    handles.coefficients=present_coefficients(output.coeff);   
    handles.Max_Dev=output.Max_Dev;
    handles.zfitted=output.zfitted;
    set(handles.results,'String',mat2str(handles.coefficients,5)); 
    set(handles.plot,'Enable','on');
    set(handles.accept,'Enable','on');
    set(handles.exportdxf,'Enable','on');
    guidata(h, handles);
    update_dev_graph(handles);
else
    if handles.maxdevvalue==[]
        handles.maxdevvalue=0
    else
    end
    if (get(handles.allowselected,'Value'))==1
            output=feval(@optim_kernel,handles.dataset,handles.series,handles.maxdevvalue,handles.bcx_0,handles.bcx_end,handles.bcy_0,handles.bcy_end);
    else
            output1=feval(@optim_kernel,handles.dataset,'SIN*SIN',handles.maxdevvalue,handles.bcx_0,handles.bcx_end,handles.bcy_0,handles.bcy_end);
            output2=feval(@optim_kernel,handles.dataset,'(COS^2)*(COS^2)',handles.maxdevvalue,handles.bcx_0,handles.bcx_end,handles.bcy_0,handles.bcy_end);
        if output1.Max_Dev>output2.Max_Dev
                output=output2;
                set(handles.seriesselector,'Value',2);
        else
                output=output1;
                set(handles.seriesselector,'Value',1);
        end
    end
  handles.coefficients=present_coefficients(output.coeff);   
  handles.Max_Dev=output.Max_Dev;
  handles.order=output.optim_order; guidata(h, handles); 
  handles.zfitted=output.zfitted;
  set(handles.results,'String',mat2str(handles.coefficients,5)); 
  set(handles.plot,'Enable','on');
  set(handles.accept,'Enable','on');
  set(handles.exportdxf,'Enable','on');
  update_dev_graph(handles); 
  
  handles.order=3; guidata(h, handles);
  
end
    

% --------------------------------------------------------------------
function varargout = dataselector_Callback(h, eventdata, handles, varargin)
dataset_name=get_var_names(handles);
handles.dataset_name=dataset_name;
handles.dataset=evalin('base',dataset_name);  guidata(h, handles);
dataset_size=size(handles.dataset);
if (dataset_size(2))<3
    set(handles.text2,'Enable','off');
    set(handles.results,'String','target dataset must contain at least 3 columns...');
    handles.dataset=[]; guidata(h, handles);
else
    set(handles.results,'String','target dataset OK...');
    set(handles.text2,'Enable','on');
    set(handles.text4,'Enable','on');
    set(handles.plot,'Enable','off');
    set(handles.accept,'Enable','off');
    set(handles.exportdxf,'Enable','off');
end


% --------------------------------------------------------------------
function varargout = seriesselector_Callback(h, eventdata, handles, varargin)
handles.series=get_series_selection(handles);
guidata(h, handles);


% --------------------------------------------------------------------
function varargout = bcx0f_Callback(h, eventdata, handles, varargin)
off =[handles.bcx0s];
value=abs(get(off,'Value')-1);
set(off,'Value',value);
guidata(h, handles);


% --------------------------------------------------------------------
function varargout = bcx0s_Callback(h, eventdata, handles, varargin)
off =[handles.bcx0f];
value=abs(get(off,'Value')-1);
set(off,'Value',value);
guidata(h, handles);


% --------------------------------------------------------------------
function varargout = bcxmf_Callback(h, eventdata, handles, varargin)
off =[handles.bcxms];
value=abs(get(off,'Value')-1);
set(off,'Value',value);
guidata(h, handles);


% --------------------------------------------------------------------
function varargout = bcxms_Callback(h, eventdata, handles, varargin)
off =[handles.bcxmf];
value=abs(get(off,'Value')-1);
set(off,'Value',value);
guidata(h, handles);


% --------------------------------------------------------------------
function varargout = bcy0f_Callback(h, eventdata, handles, varargin)
off =[handles.bcy0s];
value=abs(get(off,'Value')-1);
set(off,'Value',value);
guidata(h, handles);


% --------------------------------------------------------------------
function varargout = bcy0s_Callback(h, eventdata, handles, varargin)
off =[handles.bcy0f];
value=abs(get(off,'Value')-1);
set(off,'Value',value);
guidata(h, handles);


% --------------------------------------------------------------------
function varargout = bcymf_Callback(h, eventdata, handles, varargin)
off =[handles.bcyms];
value=abs(get(off,'Value')-1);
set(off,'Value',value);
guidata(h, handles);


% --------------------------------------------------------------------
function varargout = bcyms_Callback(h, eventdata, handles, varargin)
off =[handles.bcymf];
value=abs(get(off,'Value')-1);
set(off,'Value',value);
guidata(h, handles);


% --------------------------------------------------------------------
function varargout = optimization_Callback(h, eventdata, handles, varargin)
state=get(handles.optimization,'Value');
if state==0
    set(handles.text3,'Enable','off');
    set(handles.maxdev,'Enable','off');
    set(handles.maxdevpc,'Enable','off');
    set(handles.allowselected,'Enable','off');
    set(handles.allowany,'Enable','off'); guidata(h, handles);
else
    set(handles.text3,'Enable','on');
    set(handles.maxdev,'Enable','on');
    set(handles.maxdevpc,'Enable','on');
    set(handles.allowselected,'Enable','on');
    set(handles.allowany,'Enable','on'); guidata(h, handles);
end


% --------------------------------------------------------------------
function varargout = maxdev_Callback(h, eventdata, handles, varargin)
valuetxt=get(handles.maxdev,'String');
value=str2num(valuetxt);
if value==[]
    set(handles.results,'String','please input a numerical value...');
    set(handles.maxdev,'String','--')
    set(handles.maxdevpc,'String','--')
    handles.maxdevvalue=Inf; guidata(h, handles);
else
    set(handles.results,'String','parameter accepted...');
    handles.maxdevvalue=value; guidata(h, handles);
    maxdeflection=max(abs(handles.dataset(:,3)));
    valuepc=100*(value/maxdeflection);
    set(handles.maxdevpc,'String',num2str(valuepc));
    
end   
    

% --------------------------------------------------------------------
function varargout = maxdevpc_Callback(h, eventdata, handles, varargin)
valuetxt=get(handles.maxdevpc,'String');
value=str2num(valuetxt);
if value==[]
    set(handles.results,'String','please input a numerical value...');
    set(handles.maxdevpc,'String','--')
    set(handles.maxdev,'String','--')
else
    set(handles.results,'String','parameter accepted...');
    maxdeflection=max(abs(handles.dataset(:,3)));
    handles.maxdevvalue=0.01*value*maxdeflection; guidata(h, handles);
    set(handles.maxdev,'String',num2str(handles.maxdevvalue));
end   


% --------------------------------------------------------------------
function varargout = allowselected_Callback(h, eventdata, handles, varargin)
off =[handles.allowany];
value=abs(get(off,'Value')-1);
set(off,'Value',value);
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = allowany_Callback(h, eventdata, handles, varargin)
off =[handles.allowselected];
value=abs(get(off,'Value')-1);
set(off,'Value',value);
guidata(h, handles);

% --------------------------------------------------------------------
function varargout = orderselector_Callback(h, eventdata, handles, varargin)
ordertxt=get_order_selection(handles);
handles.order=str2num(ordertxt);
guidata(h, handles);



% --------------------------------------------------------------------
function varargout = exportdxf_Callback(h, eventdata, handles, varargin)
x=handles.dataset(:,1); y=handles.dataset(:,2);
xline=linspace(min(x),max(x),50);
yline=linspace(min(y),max(y),50*floor(abs((max(x)-min(x))/(max(y)-min(y)))));
[X,Y]=meshgrid(xline,yline);
Zfitted=griddata(x,y,handles.zfitted,X,Y,'v4');
feval(@writedxf,handles.dataset_name,X,Y,Zfitted);