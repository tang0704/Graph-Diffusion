function varargout = test_sub_seeds(varargin)
% TEST_SUB_SEEDS MATLAB code for test_sub_seeds.fig
%      TEST_SUB_SEEDS, by itself, creates a new TEST_SUB_SEEDS or raises the existing
%      singleton*.
%
%      H = TEST_SUB_SEEDS returns the handle to a new TEST_SUB_SEEDS or the handle to
%      the existing singleton*.
%
%      TEST_SUB_SEEDS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEST_SUB_SEEDS.M with the given input arguments.
%
%      TEST_SUB_SEEDS('Property','Value',...) creates a new TEST_SUB_SEEDS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before test_sub_seeds_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to test_sub_seeds_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help test_sub_seeds

% Last Modified by GUIDE v2.5 31-Oct-2013 16:44:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @test_sub_seeds_OpeningFcn, ...
                   'gui_OutputFcn',  @test_sub_seeds_OutputFcn, ...
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
addpath(genpath('.'));





% --- Executes just before test_sub_seeds is made visible.
function test_sub_seeds_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to test_sub_seeds (see VARARGIN)

% Choose default command line output for test_sub_seeds
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes test_sub_seeds wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = test_sub_seeds_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in priortype.
function priortype_Callback(hObject, eventdata, handles)
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
handles.prior_type = 'ground';
guidata(hObject,handles);

% --- Executes on button press in open.
function open_Callback(hObject, eventdata, handles)

handles.database = 'MSRA1000'; % MSRA1000, Berkeley300
handles.database_out = 'MSRA1000';
sp_num_max = 200;
graphOpts.featMode=1; 
if graphOpts.featMode == 3
    graphOpts.feat_dist_opt =  'SDMVC';
end

%%
[imdir, handles.spdir, handles.saldir,handles.gdir] = gene_dir(handles.database,handles.database_out);

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

%% compute the feature (mean color in lab color space)        
[sp_fea, rgb_fea]  = gene_feature(handles.input_im, handles.superpixels, handles.sp_center, handles.sp_npix, graphOpts);        
min_sp_fea = min(sp_fea(:)); min_rgb_fea = min(rgb_fea(:));
handles.sp_fea=(sp_fea-min_sp_fea)/(max(sp_fea(:))-min_sp_fea);
handles.rgb_fea=(rgb_fea-min_rgb_fea)/(max(rgb_fea(:))-min_rgb_fea);
%%
handles.init_seeds = [];
handles.candi_seeds = [];
% handles.prior_type = 'ground';
% set(handles.priortype,'String',handles.prior_type);
% contents = cellstr(get(hObject,'String'));
% handles.priortype = contents{get(hObject,'Value')};
guidata(hObject,handles)

% --- Executes on button press in initSeeds.
function initSeeds_Callback(hObject, eventdata, handles)
[x, y] = getpts(handles.axes1);
for i=1:length(x)
    tmp = handles.superpixels(round(y(i)),round(x(i)));
    handles.init_seeds = update_seeds(tmp, handles.init_seeds);
end
draw_seeds(handles,handles.axes1);
guidata(hObject,handles);

function seeds_out = update_seeds(sp_no, seeds_in)
seeds_out = seeds_in;
rr = find( seeds_out == sp_no);
if isempty(rr)
    seeds_out = [seeds_out; sp_no];
else
    seeds_out(rr) = [];     
end

function draw_seeds(handles, axesi, i)

initSeedColor = zeros(length(handles.init_seeds),3);
initSeedColor(:,1) = 255;

if nargin == 3
    candi_seeds = handles.candi_seeds(i);
else
    candi_seeds = handles.candi_seeds;
end

candiSeedColor = zeros(length(candi_seeds),3);
candiSeedColor(:,2) = 255;
seeds = [handles.init_seeds; candi_seeds];
seedcolor = [initSeedColor; candiSeedColor];

sp_im = draw_seed(handles.input_im, handles.superpixels, seeds, seedcolor);    
axes(axesi);
imshow(sp_im); 


% --- Executes on button press in candiSeeds.
function candiSeeds_Callback(hObject, eventdata, handles)
[x, y] = getpts(handles.axes1);
for i=1:length(x)
    tmp = handles.superpixels(round(y(i)),round(x(i)));
    handles.candi_seeds = update_seeds(tmp, handles.candi_seeds);
end

draw_seeds(handles,handles.axes1);
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
graphOpts.second_graph_method = 'feature';
%% parameters group for selecting seeds
submodularOpts.ground_conductance_coeff_foreground = 0.01; 
submodularOpts.ground_conductance_coeff_background = 1;  

submodularOpts.alpha = 0.01; 
submodularOpts.Walpha = 0;   %0 : (D - W +alpha*I)F = alpha*Y  1:  (D - 1/(1+alpha)*W)F = alpha/(1+alpha)*D*Y. Both are same!!!

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
        if strcmp(handles.database,'MSRA1000')
            gname=[handles.gdir handles.imname(1:end-4) '_G.bmp'];
        else
            gname=[handles.gdir handles.imname(1:end-4) '.png'];
        end
          
            
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
        if strcmp(handles.database,'MSRA1000')
            gname=[handles.gdir handles.imname(1:end-4) '_G.bmp'];
        else
            gname=[handles.gdir handles.imname(1:end-4) '.png'];
        end
        
        [neg_label_inds, EBW_sp, CBW_sp,BW_sp] =  gene_convex_curve(spAdjcMat, superpixels,handles.sp_center,...
                                                  handles.input_im, priorOpts,1,handles.saldir,handles.imname);
        close; close;close;
       %% compute 3 priors: location, color & background
        cueImgCnt = cue_by_img_center(handles.sp_center, superpixels);%sp_center, sp_npix            
        tmp = min(cueImgCnt); cueImgCnt=(cueImgCnt-tmp)/(max(cueImgCnt(:))-tmp);
        [c1, cue_by_color] = cue_by_img_color(train_data_dir, superpixels,...
                             handles.sp_center, handles.sp_inds,rgb_fea,[]); 
        tmp = min(cue_by_color); cue_by_color=(cue_by_color-tmp)/(max(cue_by_color(:))-tmp);

        bd=unique([superpixels(1,:),superpixels(m,:),superpixels(:,1)',superpixels(:,n)']); 
        bd  = setdiff(bd',neg_label_inds);
        bd = [neg_label_inds;bd]; 
         
        aff_mat = gene_weight( spAdjcMat, sp_fea, theta, graphOpts,bd,[],2);% affinity matrix,   

        [prior_source] = gene_prior_source(aff_mat, neg_label_inds, priorOpts.doBinary, graphOpts);
         min_c = min(prior_source);
         prior_source = (prior_source - min_c)./(max(prior_source) - min_c);
   
        prior_source_with_cue = prior_source.*cueImgCnt.*cue_by_color;

        minpc = min(prior_source_with_cue(:));
        prior_source_with_cue=(prior_source_with_cue-minpc)/(max(prior_source_with_cue)-minpc);
end
axes(handles.axes13);
mcolor = zeros(length(neg_label_inds),3);
mcolor(:,3) = 144;
% [label_im_neg] = draw_seed(handles.input_im, superpixels, neg_label_inds, mcolor);
% poutname=[handles.saldir handles.imname(1:end-4) '_bseed.png'];
% imwrite(label_im_neg,poutname);
poutname=[handles.saldir handles.imname(1:end-4) '_foreground_convex.png'];
foreconvex = imread(poutname);
imshow(foreconvex);
% imshow(label_im_neg);

axes(handles.axes14);
[prior_wiht_cue_map] = saliency_sp2im(prior_source_with_cue, ...
            show_para.sp_inds, show_para.sp_num, show_para.m, show_para.n, show_para.w);    
poutname=[handles.saldir handles.imname(1:end-4) '_p_cue.png'];
imwrite(prior_wiht_cue_map,poutname);
imshow(prior_wiht_cue_map);

%% perform learning-based les / select seeds!  

% bd=unique([superpixels(1,:),superpixels(m,:),superpixels(:,1)',superpixels(:,n)']);
% bd  = setdiff(bd',neg_label_inds);
% bd = [neg_label_inds;bd]; 
% aff_mat = gene_weight( spAdjcMat,[sp_fea  repmat(prior_source_with_cue,1,3)], theta, graphOpts, bd,[],2);
switch graphOpts.second_graph_method,
    case 'feature'
        aff_mat = gene_weight( spAdjcMat, sp_fea, theta, graphOpts, bd,[],2);
    case 'multi'
       aff_mat = gene_weight( spAdjcMat, [sp_fea  repmat(prior_source,1,3)], theta, graphOpts, bd,[],2);
end
% aff_mat = gene_weight( spAdjcMat, sp_fea, theta, graphOpts, neg_label_inds,[],2);

for i=1:min(5, length(handles.candi_seeds))
    pos_label_inds = [handles.init_seeds;handles.candi_seeds(i)];
    pos_label_ranks = prior_source_with_cue(pos_label_inds);   
    % sp_saliency = solve_les_with_dirichlet_no_prior(aff_mat, pos_label_ranks, pos_label_inds, graphOpts);
    sp_saliency = solve_les_with_dirichlet_v2(aff_mat, pos_label_inds, pos_label_ranks,neg_label_inds, prior_source_with_cue, submodularOpts);
    pos_label_inds_choosen = pos_label_inds;
%     tmp = min(sp_saliency);
%     show_sp_saliency=(sp_saliency-tmp)/(max(sp_saliency)-tmp);
%     L(i)=sum(show_sp_saliency);
        L(i)=sum(sp_saliency);
%     LW(i) = L(i) - sum(1./(pos_label_ranks.^3));  

    %% assign the saliency value to each pixel
    [label_im] = draw_seed(handles.input_im, handles.superpixels,  pos_label_inds_choosen);
    soutname = sprintf('%s-%d_seed.png', [handles.saldir handles.imname(1:end-4)], i);
    imwrite(label_im,soutname);
    [im_saliency] = saliency_sp2im(sp_saliency, handles.sp_inds, handles.sp_num, m, n, w,0);

    % post process: just smoothing!
    myfilter = fspecial('gaussian',[5 5], 1);
    im_saliency = imfilter(im_saliency, myfilter, 'replicate');

    %% output
    switch i
        case 1
            axis1 = handles.axes3;
            axis2 = handles.axes8;
        case 2
            axis1 = handles.axes4;
            axis2 = handles.axes9;
        case 3
            axis1 = handles.axes5;
            axis2 = handles.axes10;
        case 4
            axis1 = handles.axes6;
            axis2 = handles.axes11;            
        case 5
            axis1 = handles.axes7;
            axis2 = handles.axes12;            
    end
    draw_seeds(handles,axis1,i);
    axes(axis2);    imshow(im_saliency);

    outname = sprintf('%s-%d.bmp', [handles.saldir handles.imname(1:end-4)], i);    
    imwrite(im_saliency,outname);
end

%% submodular
pos_label_inds = find(BW_sp);
% pos_label_inds = gene_pos_retrive_pts(prior_source_with_cue, BW_sp, fseedOpts);     
seed_num = numel(pos_label_inds);
submodularOpts.gamma ='prior';
%  aff_mat = gene_weight( spAdjcMat, [sp_fea  repmat(prior_source,1,3)], theta, graphOpts, bd,[],2);
% 121,,74,superpixels(195,154),superpixels(222,150),71

% pos_label_inds = [ 121,71,superpixels(195,154),74];    
% pos_label_ranks = prior_source_with_cue(pos_label_inds);
% sp_saliency = solve_les_with_dirichlet_v2(aff_mat, pos_label_inds, pos_label_ranks,neg_label_inds, prior_source_with_cue, submodularOpts);
% pos_label_inds_choosen = pos_label_inds;

[sp_saliency, pos_label_inds_choosen] = solve_submodular_naive(...
    aff_mat, seed_num, pos_label_inds, neg_label_inds, prior_source_with_cue, ...
    submodularOpts,handles.input_im, superpixels,show_para,0,handles.saldir,handles.imname); 

% tmp = min(sp_saliency);
% show_sp_saliency=(sp_saliency-tmp)/(max(sp_saliency)-tmp);

L(i+1)=sum(sp_saliency);
% LW(i+1) = L(i+1) - sum(1./(prior_source_with_cue(handles.candi_seeds(i)).^3));
% LW(i+1) = L(i+1) - sum(1./(prior_source_with_cue(pos_label_inds_choosen).^2));    
   
%% assign the saliency value to each pixel
[label_im] = draw_seed(handles.input_im, handles.superpixels,  pos_label_inds_choosen);
soutname = sprintf('%s.png', [handles.saldir handles.imname(1:end-4) '_our_seed']);
imwrite(label_im,soutname);
axes(handles.axes15);
imshow(label_im); 

[im_saliency] = saliency_sp2im(sp_saliency, handles.sp_inds, handles.sp_num, m, n, w,0);
% post process: just smoothing!
myfilter = fspecial('gaussian',[5 5], 1);
im_saliency = imfilter(im_saliency, myfilter, 'replicate');
outname = sprintf('%s.bmp', [handles.saldir handles.imname(1:end-4) '_our_sal'] );    
imwrite(im_saliency,outname);
axes(handles.axes16);
imshow(im_saliency); 
%%
axes(handles.axes2); hold off;
X=1:length(L);
% plot(X, L, 'r--', 'linewidth', 1 ); hold on;
bar(X,L,'b');

legend('L');
figure;bar(X,L,'b');
saveas(gcf,[handles.saldir handles.imname(1:end-4),'_hist.png']);
L
close;
filename = [handles.saldir handles.imname(1:end-4),'_L.txt'];
fopen(filename,'w');
dlmwrite(filename, L);

% saveas(gcf,[handles.saldir handles.imname(1:end-4),'_hist.png']);



