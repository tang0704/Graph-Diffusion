function varargout = test_pde_seeds(varargin)
% TEST_PDE_SEEDS MATLAB code for test_pde_seeds.fig
%      TEST_PDE_SEEDS, by itself, creates a new TEST_PDE_SEEDS or raises the existing
%      singleton*.
%
%      H = TEST_PDE_SEEDS returns the handle to a new TEST_PDE_SEEDS or the handle to
%      the existing singleton*.
%
%      TEST_PDE_SEEDS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEST_PDE_SEEDS.M with the given input arguments.
%
%      TEST_PDE_SEEDS('Property','Value',...) creates a new TEST_PDE_SEEDS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before test_pde_seeds_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to test_pde_seeds_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help test_pde_seeds

% Last Modified by GUIDE v2.5 31-Oct-2013 13:52:44

% Begin initialization code - DO NOT EDIT
addpath(genpath('.'));
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @test_pde_seeds_OpeningFcn, ...
                   'gui_OutputFcn',  @test_pde_seeds_OutputFcn, ...
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


% --- Executes just before test_pde_seeds is made visible.
function test_pde_seeds_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to test_pde_seeds (see VARARGIN)

% Choose default command line output for test_pde_seeds
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes test_pde_seeds wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = test_pde_seeds_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in open.
function open_Callback(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

database = 'MSRA1000'; % MSRA1000, Berkeley300
sp_num_max = 200;
graphOpts.featMode=1; % bad results for 3 
if graphOpts.featMode == 3
    graphOpts.feat_dist_opt =  'SDMVC';
end

%%
[imdir, handles.spdir, handles.saldir,handles.gdir] = gene_dir(database);

[handles.imname, handles.imdir] = uigetfile({'*.bmp', 'Image Files'}, 'Pick a file', imdir);

[handles.input_im] = imread([handles.imdir handles.imname]);
[m,n,k] = size(handles.input_im); 
w = [m,n,1,m,1,n];

%%
[handles.superpixels, handles.spAdjcMat, handles.sp_inds, handles.sp_center, handles.sp_npix] = ... 
    gene_superpixel(handles.imdir, handles.imname, sp_num_max, handles.spdir, m, n);
handles.graphOpts.npix = handles.sp_npix;
handles.sp_num = size(handles.spAdjcMat,1);
    
axes(handles.axes1);
sp_im = draw_seed(handles.input_im, handles.superpixels, [],[255,0,0]);    
imshow(sp_im);  

handles.sp_seeds = [0;0;0];
handles.seedcolor = [255,0,0; 0, 255, 0; 0,0,255];
%% compute the feature (mean color in lab color space)        
[sp_fea, rgb_fea]  = gene_feature(handles.input_im, handles.superpixels, handles.sp_center, handles.sp_npix, graphOpts);        
min_sp_fea = min(sp_fea(:)); min_rgb_fea = min(rgb_fea(:));
handles.sp_fea=(sp_fea-min_sp_fea)/(max(sp_fea(:))-min_sp_fea);
handles.rgb_fea=(rgb_fea-min_rgb_fea)/(max(rgb_fea(:))-min_rgb_fea);

guidata(hObject,handles)

% --- Executes on button press in seed1.
function seed1_Callback(hObject, eventdata, handles)
[x, y] = getpts(handles.axes1);
tmp = handles.superpixels(round(y),round(x));
handles = update_seeds(1, tmp, handles);
guidata(hObject,handles);

function handles = update_seeds(i, sp_no, handles)
rr = find( handles.sp_seeds == sp_no);
if isempty(rr)
    handles.sp_seeds(i) = sp_no;
elseif rr == i
    handles.sp_seeds(i) = 0;       
else
    handles.sp_seeds(i) = sp_no;
    handles.sp_seeds(rr) = 0;    
end
sp_im = draw_seed(handles.input_im, handles.superpixels, handles.sp_seeds, handles.seedcolor);    
axes(handles.axes1);
imshow(sp_im); 

% --- Executes on button press in seed2.
function seed2_Callback(hObject, eventdata, handles)
[x, y] = getpts(handles.axes1);
tmp = handles.superpixels(round(y),round(x));
handles = update_seeds(2, tmp, handles);
guidata(hObject,handles);

% --- Executes on button press in seed3.
function seed3_Callback(hObject, eventdata, handles)
[x, y] = getpts(handles.axes1);
tmp = handles.superpixels(round(y),round(x));
handles = update_seeds(3, tmp, handles);
guidata(hObject,handles);
   
% --- Executes on button press in compute.
function compute_Callback(hObject, eventdata, handles)
theta = 10; % control the edge weight
train_data_dir = 'TrainedData';
%% back & front seed parameters, 用缩进表示稳定的参数
    priorOpts.enlarge_convex_ratio = 0.1;   
    priorOpts.bseed_thresh = 26;     
    priorOpts.bseed_percentage =1;        
        priorOpts.doBinary = 0; %0: do not binary the prior 1: binary the prior  (0 is better)

    fseedOpts.fseed_percentage = 1;   

%% graph parameters
graphOpts.ground_conductance_coeff = 0.1; 

%% parameters group for selecting seeds
submodularOpts.ground_conductance_coeff_foreground = 0.01; 
submodularOpts.ground_conductance_coeff_background = 1;  

submodularOpts.alpha = 0.01; 
submodularOpts.Walpha = 1;   %0 : (D - W +alpha*I)F = alpha*Y  1:  (D - 1/(1+alpha)*W)F = alpha/(1+alpha)*D*Y. Both are same!!!

submodularOpts.seedtype = 'initial'; % initial: prior map as the seed rank;  ones: [1...1,0]  as the seed rank; current: previous saliency as the seed rank
submodularOpts.ytype = 'initial'; % initial: prior map as  Y  current: previous saliency map  as  Y

submodularOpts.gamma = 'prior'; %  control the seed number, 'prior' or a number

%%
sp_fea = handles.sp_fea;
rgb_fea = handles.rgb_fea;
superpixels = handles.superpixels;
spAdjcMat = handles.spAdjcMat;
[m,n,k] = size(handles.input_im); 
w = [m,n,1,m,1,n];

show_para = [];
show_para.sp_inds = handles.sp_inds;    show_para.sp_num = handles.sp_num;
show_para.m = m;    show_para.n = n;    show_para.w = w; 
    
switch handles.prior_type
    case 'ground'        
        gname=[handles.gdir handles.imname(1:end-4) '_G.bmp'];
        prior_map = imread(gname);
        prior_source_with_cue = img2superpixel(prior_map,superpixels)./255;
        
        neg_label_inds = find(prior_source_with_cue==0);
        
        BW_sp = (prior_source_with_cue==1);
        EBW_sp = BW_sp;
        tmp = find(BW_sp == 1);
        for i = tmp'
            EBW_sp(handles.spAdjcMat(i,:)>0) = 1;
        end
        
        CBW_sp = BW_sp;
        tmp = find(BW_sp == 0);
        for i = tmp'
            CBW_sp(handles.spAdjcMat(i,:)>0) = 0;
        end
    otherwise
        % generate background seeds    
        [neg_label_inds, EBW_sp, CBW_sp, BW_sp] = gene_neg_retrive_pts(...
                spAdjcMat, superpixels, handles.input_im, priorOpts, ...
                0, 0, handles.saldir,handles.imname);

       %% compute 3 priors: location, color & background
        cueImgCnt = cue_by_img_center(handles.sp_center, superpixels);%sp_center, sp_npix            
        tmp = min(cueImgCnt); cueImgCnt=(cueImgCnt-tmp)/(max(cueImgCnt(:))-tmp);
        %             cueImgCnt = ones(sp_num,1);
        [c1, cue_by_color] = cue_by_img_color(train_data_dir, superpixels,...
                             handles.sp_center, handles.sp_inds,rgb_fea,[]); 
        tmp = min(cue_by_color); cue_by_color=(cue_by_color-tmp)/(max(cue_by_color(:))-tmp);

        bd=unique([superpixels(1,:),superpixels(m,:),superpixels(:,1)',superpixels(:,n)']); 
        bd  = setdiff(bd',neg_label_inds);
        bd = [neg_label_inds;bd]; 
        %     bd = []; 
        aff_mat = gene_weight( spAdjcMat, sp_fea, theta, graphOpts,bd,[],2);% affinity matrix,   
        [prior_source] = gene_prior_source(aff_mat, neg_label_inds, priorOpts.doBinary, graphOpts);

        %     
        %             tmp = cueImgCnt.*cue_by_color;
        %             prior_source_with_cue = priorOpts.alpha*prior_source + (1-priorOpts.alpha)*(tmp);
        prior_source_with_cue = prior_source.*cueImgCnt.*cue_by_color;

        minpc = min(prior_source_with_cue(:));
        prior_source_with_cue=(prior_source_with_cue-minpc)/(max(prior_source_with_cue)-minpc);
end
axes(handles.axes2);
mcolor = zeros(length(neg_label_inds),3);
mcolor(:,3) = 144;
[label_im_neg] = draw_seed(handles.input_im, superpixels, neg_label_inds, mcolor);
poutname=[handles.saldir handles.imname(1:end-4) '_bseed.png'];
imwrite(label_im_neg,poutname);
imshow(label_im_neg);

axes(handles.axes3);
[prior_wiht_cue_map] = saliency_sp2im(prior_source_with_cue, ...
            show_para.sp_inds, show_para.sp_num, show_para.m, show_para.n, show_para.w);    
poutname=[handles.saldir handles.imname(1:end-4) '_p_cue.png'];
imwrite(prior_wiht_cue_map,poutname);
imshow(prior_wiht_cue_map);

%% perform learning-based les / select seeds!  算完初始显著值后重新定义图
pos_label_inds = handles.sp_seeds(handles.sp_seeds>0);
pos_label_ranks = prior_source_with_cue(pos_label_inds);   
pos_label_ranks(pos_label_ranks<1) = 0.5;

bd=unique([superpixels(1,:),superpixels(m,:),superpixels(:,1)',superpixels(:,n)']);
bd  = setdiff(bd',neg_label_inds);
bd = [neg_label_inds;bd]; 
aff_mat = gene_weight( spAdjcMat,[sp_fea  repmat(prior_source_with_cue,1,3)], theta, graphOpts, bd,[],2);
% aff_mat = gene_weight( spAdjcMat, sp_fea, theta, graphOpts, bd,[],2);

% sp_saliency = solve_les_with_dirichlet_no_prior(aff_mat, pos_label_ranks, pos_label_inds, graphOpts);
sp_saliency = solve_les_with_dirichlet_v2(aff_mat, pos_label_inds, pos_label_ranks,neg_label_inds, prior_source_with_cue, submodularOpts);
pos_label_inds_choosen = pos_label_inds;

tmp = min(sp_saliency);
show_sp_saliency=(sp_saliency-tmp)/(max(sp_saliency)-tmp);

L=sum(show_sp_saliency)
tmp = sprintf('%f', L);
set(handles.textL,'String',tmp);

seed_num = numel(pos_label_inds);
submodularOpts.gamma = seed_num*0.2;
LW = L - submodularOpts.gamma *sum(1./pos_label_ranks.^2);    
tmp = sprintf('%f', LW);
set(handles.textLW,'String',tmp);

%% assign the saliency value to each pixel
[label_im] = draw_seed(handles.input_im, handles.superpixels,  pos_label_inds_choosen, [255,0,0]);
soutname=[handles.saldir handles.imname(1:end-4) '.png'];
imwrite(label_im,soutname);
[im_saliency] = saliency_sp2im(sp_saliency, handles.sp_inds, handles.sp_num, m, n, w);

% post process: just smoothing!
myfilter = fspecial('gaussian',[5 5], 1);
im_saliency = imfilter(im_saliency, myfilter, 'replicate');

%% output
axes(handles.axes4);
imshow(im_saliency);

outname=[handles.saldir handles.imname(1:end-4) '.bmp'];
imwrite(im_saliency,outname);


% --- Executes on selection change in priortype.
function priortype_Callback(hObject, eventdata, handles)
% hObject    handle to priortype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns priortype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from priortype
contents = cellstr(get(hObject,'String'));
handles.prior_type = contents{get(hObject,'Value')};
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function priortype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to priortype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
