function varargout = photoshop(varargin)
% PHOTOSHOP MATLAB code for photoshop.fig
%      PHOTOSHOP, by itself, creates a new PHOTOSHOP or raises the existing
%      singleton*.
%
%      H = PHOTOSHOP returns the handle to a new PHOTOSHOP or the handle to
%      the existing singleton*.
%
%      PHOTOSHOP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHOTOSHOP.M with the given input arguments.
%
%      PHOTOSHOP('Property','Value',...) creates a new PHOTOSHOP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before photoshop_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to photoshop_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help photoshop

% Last Modified by GUIDE v2.5 13-Dec-2017 21:25:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @photoshop_OpeningFcn, ...
                   'gui_OutputFcn',  @photoshop_OutputFcn, ...
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


% --- Executes just before photoshop is made visible.
function photoshop_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to photoshop (see VARARGIN)

% Choose default command line output for photoshop
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes photoshop wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = photoshop_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadImage.
function loadImage_Callback(hObject, eventdata, handles)
% hObject    handle to loadImage (see GCBO) 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image_in threshold image_out;
[filename , pathname] = uigetfile({'*.jpg';'*.png'},'File Selector');
image_in = imread(strcat(pathname, filename));
aspectRatio= size(image_in,1)/size(image_in,2);
if aspectRatio > 1
    x=512 - size(image_in,1);
    image_in=imresize(image_in,[size(image_in,1)+x size(image_in,2)+x/aspectRatio]);
else
    x=512 - size(image_in,2);
    image_in=imresize(image_in,[size(image_in,1)+x*aspectRatio size(image_in,2)+x]);
end
image_in=RGBtoGray(image_in);
axes(handles.axes1);
imshow(image_in);
image_out = image_in;
threshold = getGlobalThreshold(image_in);

set(handles.LaplacianSharpen,'Enable','on');%Turns button two on when image is loaded.
set(handles.histogramEqualization,'Enable','on');
set(handles.globalThreshold,'Enable','on');
set(handles.gaussianBlur,'Enable','on');
set(handles.prewittEdge,'Enable','on');
set(handles.sobelEdge,'Enable','on');
set(handles.zoom,'Enable','on');
set(handles.dithering,'Enable','on');
set(handles.saveImage,'Enable','on');
set(handles.rotateImageDegrees,'Enable','on');


function Y = RGBtoGray(image)
r=image(:,:,1);g=image(:,:,2);b=image(:,:,3);
ITU=[0.299,0.587,0.114];
Y=round(r*ITU(1)+g*ITU(2)+b*ITU(3));

function t = getGlobalThreshold(image)
a=5;
t= floor(rand * 256);
ti= 300;

while (abs(t-ti) >= a)    
    ti=t;
    above_t= image>t;

    sumG1= sum(image(above_t));
    sizeG1= size(image(above_t));
    
    sumG2= sum(image(~above_t));
    sizeG2= size(image(~above_t));
    
    t=((sumG1/sizeG1(1))+(sumG2/sizeG2(1)))/2;
end

% --- Executes on button press in LaplacianSharpen.
function LaplacianSharpen_Callback(hObject, eventdata, handles)
%RGB to Grayscale
% hObject    handle to LaplacianSharpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image_in image_out;
LP=[0,-1,0;
   -1,5,-1;
    0,-1,0];
image_out = filter3by3(LP, image_in);
axes(handles.axes2);
imshow(image_out);

function res= filter3by3(filter,image)
res = image;
image = double(image);
sz = size(image);
for j= 2:sz(1)-1
    for k= 2:sz(2)-1
        f=sum(sum(image(j-1:j+1,k-1:k+1) .* filter));
        if f >255
            f = 255;
        elseif f < 0
            f = 0;
        end
        res(j,k)=f;
    end
end

% --- Executes on button press in histogramEqualization.
function histogramEqualization_Callback(hObject, eventdata, handles)
% hObject    handle to histogramEqualization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image_in image_out;
image_out=image_in;
copy_image = image_in;
pixelFrequency = double(zeros(1,256));
s = double(zeros(1,256));
sum = double(0);
sz = size(copy_image);

for i= 1:sz(1)
    for j= 1:sz(2)
       i=uint32(i);
        j=uint32(j);
       pixelFrequency(copy_image(i,j)+1) =  pixelFrequency(copy_image(i,j)+1) + 1;
    end
end

for i= 1:256
    i=uint32(i);
       pixelFrequency(i) =  pixelFrequency(i)/(sz(1)*sz(2));
       sum = sum + pixelFrequency(i);
       s(i) = round(255*sum);
end

for i= 1:sz(1)
    for j= 1:sz(2)
       i=uint32(i);
       j=uint32(j);
       image_out(i,j)= s(copy_image(i,j)+1);
    end
end

axes(handles.axes2);
imshow(image_out);

% --- Executes on button press in globalThreshold.
function globalThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to globalThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image_in threshold image_out;
image_out=image_in;
i=image_in>threshold;
image_out(i) = 255;
image_out(~i) = 0;
axes(handles.axes2);
imshow(image_out);

% --- Executes on button press in prewittedge.
function gaussianBlur_Callback(hObject, eventdata, handles)
% hObject    handle to prewittedge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image_in image_out;
edge=[1/9,1/9,1/9;
      1/9,1/9,1/9;
      1/9,1/9,1/9];
image_out = filter3by3(edge, image_in);
axes(handles.axes2);
imshow(image_out);

% --- Executes on button press in prewittEdge.
function prewittEdge_Callback(hObject, eventdata, handles)
% hObject    handle to prewittEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image_in threshold image_out;
image_out=image_in;

prewittHorizontal=[-1,-1,-1;
      0,0,0;
      1,1,1];
prewittVertical = [-1,0,1;
      -1,0,1;
      -1,0,1];
gx= filter3by3(prewittHorizontal, image_in);
gy= filter3by3(prewittVertical,image_in);

mag = abs(gx)+abs(gy);

indices = mag > threshold;

image_out(indices) = 255;
image_out(~indices)=0; 

axes(handles.axes2);
imshow(image_out);

% --- Executes on button press in sobelEdge.
function sobelEdge_Callback(hObject, eventdata, handles)
% hObject    handle to sobelEdge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image_in threshold image_out;
image_out=image_in;

sobelHorizontal=[-1,-2,-1;
      0,0,0;
      1,2,1];
sobelVertical = [-1,0,1;
      -2,0,2;
      -1,0,1];
gx= filter3by3(sobelHorizontal, image_in);
gy= filter3by3(sobelVertical,image_in);

mag = abs(gx)+abs(gy);
indices = mag > threshold;

image_out(indices) = 255;
image_out(~indices)=0; 

axes(handles.axes2);
imshow(image_out);


% --- Executes on button press in zoom.
function zoom_Callback(hObject, eventdata, handles)
% hObject    handle to zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global image_in image_out;
message=sprintf('Pick the corner points of a rectangular region to crop');
uiwait(msgbox(message));
p = ginput(2);
xmin= min(round(p(1,1)), round(p(2,1))); %Finds the Xmin coordinate
ymin = min(round(p(1,2)), round(p(2,2))); %Finds the Ymin coordinate
xmax = max(round(p(1,1)), round(p(2,1))); %Finds the Xmax coordinate
ymax = max(round(p(1,2)), round(p(2,2))); %Finds the Ymax coordinate
width = xmax-xmin+1;
height = ymax-ymin+1;
hold on;
r=rectangle('Position', [ xmin ymin width height] );
image_out=image_in(ymin:ymax,xmin:xmax,:);
axes(handles.axes2);
imshow(image_out);
waitforbuttonpress;
delete(r);


% --- Executes on button press in dithering.
function dithering_Callback(hObject, eventdata, handles)
% hObject    handle to dithering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image_in image_out;
image_out = image_in;
sz = size(image_out);
for i= 1:sz(1)-1
    for j= 2:sz(2)-1
       oldpixel=image_out(i,j);
       newpixel=255*round(oldpixel/256);
       image_out(i,j) = newpixel;
       quant_error = oldpixel-newpixel;
       image_out(i,j+1) = image_out(i,j+1) + quant_error*(7/16);
       image_out(i+1,j-1) = image_out(i+1,j-1) + quant_error*(3/16);
       image_out(i+1,j) = image_out(i+1,j) + quant_error*(5/16);
       image_out(i+1,j+1) = image_out(i+1,j+1) + quant_error*(1/16);
    end
end
axes(handles.axes2);
imshow(image_out);


% --- Executes on button press in saveImage.
function saveImage_Callback(hObject, eventdata, handles)
% hObject    handle to saveImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image_out;
imwrite(image_out,'outFile.png');
message=sprintf(strcat('Image Saved as outFile.png'));
uiwait(msgbox(message));

% --- Executes on button press in rotateImageDegrees.
function rotateImageDegrees_Callback(hObject, eventdata, handles)
% hObject    handle to rotateImageDegrees (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global image_in image_out;

prompt={'Enter a value of \theta (in degrees)'};
name = 'Theta Value';
defaultans = {'0'};
options.Interpreter = 'tex';
answer = inputdlg(prompt,name,[1 40],defaultans,options);
degree= str2double(answer);

if isnan(degree)
    image_out=image_in;
else
    degree = mod(degree,360);
    if degree == 90
        image_out = transpose(image_in);
    else

    rad=deg2rad(degree);
    copy_image = image_in;
    sz=size(copy_image);
    rowsf=ceil(sz(1)*abs(cos(rad))+sz(2)*abs(sin(rad)));
    colsf=ceil(sz(1)*abs(sin(rad))+sz(2)*abs(cos(rad)));
    image_out=uint8(zeros([rowsf colsf]));
    outsize=size(image_out);
    for i= 1: outsize(1)
        for j= 1: outsize(2)
             m1=round((i-ceil(outsize(1)/2))*cos(rad)+(j-ceil(outsize(2)/2))*sin(rad));
             m2=round(-(i-ceil(outsize(1)/2))*sin(rad)+(j-ceil(outsize(2)/2))*cos(rad));
             m1=m1+ceil(sz(1)/2);
             m2=m2+ceil(sz(2)/2);
             if(m1>=1 && m2>=1 && m1<=sz(1) && m2<=sz(2))
                image_out(i,j)= copy_image(m1,m2);
             end
        end
    end
    end
end
axes(handles.axes2);
imshow(image_out);
