% Graph diffusion for action recognition task
% modified by Tangli Chu based on
% Demo of Learning-Based Linear Elliptic System
% for Saliency Detection
% by Risheng Liu, Guangyu Zhong, JJCAO
% rsliu@dlut.edu.cn

clc;clear all;close all;
addpath(genpath('.'));
DEBUG = 1; SHOW = 1;

%% gene dir and file name
database = 'MSRA1000'; % MSRA1000, Berkeley300
database_out = 'MSRA1000';
[imdir, spdir, saldir,gdir] = gene_dir(database,database_out);
imnames=dir([imdir '*' 'bmp']);
train_data_dir = 'TrainedData';
%% parameters group 3: stable
sp_num_max = 200;    % inital superpixel number (upper bound number)
theta = 10; % control the edge weight
%% back & front seed parameters, ��������ʾ�ȶ��Ĳ���

priorOpts.priortype = 'convex'; % ground : groundtruth as prior; convex : general convex as prior;  matting: use Pan's mehod to matting the convex (todo, maybe not ok)
priorOpts.enlarge_convex_ratio = 0.1;  
priorOpts.bseed_thresh = 26;     
priorOpts.bseed_percentage =1;        
priorOpts.doBinary = 0; %0: do not binary the prior 1: binary the prior  (0 is better)

%% graph parameters
graphOpts.featMode=1; 
if graphOpts.featMode == 3
    graphOpts.feat_dist_opt =  'SDMVC';
end
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
fseedOpts.fseed_num  = 10;
%%
% Berkeley300: 42044, 105025, 216053,
% MSRA1000:
% Imgs_0978,Imgs_0977,Imgs_0964,Imgs_0937,Imgs_0932,Imgs_0924,Imgs_0915,Imgs_0907,Imgs_0905,Imgs_0894,Imgs_0891,Imgs_0886,
% Imgs_0885,Imgs_0863,Imgs_0858,Imgs_0839,Imgs_0834,Imgs_0832,Imgs_0826,Imgs_0825,Imgs_0822,Imgs_0811,Imgs_0806,Imgs_0797
% Imgs_0777,Imgs_0752,Imgs_0748,Imgs_0747,Imgs_0746,Imgs_0739
for ii=1:1
    close all;    clear w;
    %% read image and preprocessing: : just smoothing!
    imname=imnames(ii).name;
    %     if ~strcmp(imname, '8023.bmp')
    %         continue;
    %     end
    sprintf('%d: %s', ii, imname)
    [input_im] = imread([imdir, imname]);    
    [m,n,k] = size(input_im);
    w = [m,n,1,m,1,n];
    
    %% generate superpixels
    [superpixels, spAdjcMat, sp_inds, sp_center, sp_npix] = gene_superpixel(imdir, imname, sp_num_max, spdir, m, n);
    %     [superpixels, spAdjcMat, sp_inds, sp_center, sp_npix] = run_superpixel(input_im, sp_num_max) ;  
    graphOpts.npix = sp_npix;
    sp_num = size(spAdjcMat,1);
    
    %% compute the feature (mean color in lab color space)
    [sp_fea, rgb_fea]  = gene_feature(input_im, superpixels, sp_center, sp_npix, graphOpts);
    min_sp_fea = min(sp_fea(:)); min_rgb_fea = min(rgb_fea(:));
    sp_fea=(sp_fea-min_sp_fea)/(max(sp_fea(:))-min_sp_fea);
    rgb_fea=(rgb_fea-min_rgb_fea)/(max(rgb_fea(:))-min_rgb_fea);
    
    show_para = [];
    show_para.sp_inds = sp_inds;    show_para.sp_num = sp_num;
    show_para.m = m;    show_para.n = n;    show_para.w = w;
    
    %%
    switch priorOpts.priortype
        case 'ground'
            gname=[gdir imnames(ii).name(1:end-4) '_G.bmp'];
            prior_map = imread([gname]);
            prior_source_with_cue = img2superpixel(prior_map,superpixels)./255;
            neg_label_inds = [];
            CBW_sp = (prior_source_with_cue==1);
        otherwise
            % generate background seeds
            %             [neg_label_inds, EBW_sp, CBW_sp, BW_sp] = gene_neg_retrive_pts(spAdjcMat, superpixels, input_im, priorOpts, 0, SHOW,saldir,imname);
            [neg_label_inds, EBW_sp, CBW_sp, BW_sp] = gene_convex_curve(spAdjcMat, superpixels, sp_center,input_im, priorOpts,SHOW,saldir,imname);
            
            %% compute 3 priors: location, color & background
            cueImgCnt = cue_by_img_center(sp_center, superpixels);%sp_center, sp_npix
            tmp = min(cueImgCnt); cueImgCnt=(cueImgCnt-tmp)/(max(cueImgCnt)-tmp);
            [c1, cue_by_color] = cue_by_img_color(train_data_dir, superpixels,sp_center,sp_inds,rgb_fea,show_para);
            tmp = min(cue_by_color); cue_by_color=(cue_by_color-tmp)/(max(cue_by_color)-tmp);
    
            bd=unique([superpixels(1,:),superpixels(m,:),superpixels(:,1)',superpixels(:,n)']); 
            bd  = setdiff(bd',neg_label_inds);
            bd = [neg_label_inds;bd]; 
            
            aff_mat = gene_weight( spAdjcMat, sp_fea, theta, graphOpts,bd,[],2);% affinity matrix,
            [prior_source] = gene_prior_source(aff_mat, neg_label_inds, priorOpts.doBinary, graphOpts);
            tmp = min(prior_source); prior_source=(prior_source-tmp)/(max(prior_source)-tmp);
            
            prior_source_with_cue = prior_source.*cueImgCnt.*cue_by_color;
            minpc = min(prior_source_with_cue(:));
            prior_source_with_cue=(prior_source_with_cue-minpc)/(max(prior_source_with_cue)-minpc);
            
            prior_with_cue_map = saliency_sp2im(prior_source_with_cue, ...
                show_para.sp_inds, show_para.sp_num, show_para.m, ...
                show_para.n, show_para.w);
            imwrite(prior_with_cue_map,[saldir imnames(ii).name(1:end-4) '_p_cue.png']);
            if SHOW        
                prior_map = saliency_sp2im(prior_source, ...
                    show_para.sp_inds, show_para.sp_num, show_para.m, ...
                    show_para.n, show_para.w);
                imwrite(prior_map,[saldir imnames(ii).name(1:end-4) '_background_cue.png']);
                
                cueImgCnt_map = saliency_sp2im(cueImgCnt, ...
                    show_para.sp_inds, show_para.sp_num, show_para.m, ...
                    show_para.n, show_para.w);
                imwrite(cueImgCnt_map,[saldir imnames(ii).name(1:end-4) '_center_cue.png']);
                
                cue_by_color_map = saliency_sp2im(cue_by_color, ...
                    show_para.sp_inds, show_para.sp_num, show_para.m, ...
                    show_para.n, show_para.w);
                imwrite(cue_by_color_map,[saldir imnames(ii).name(1:end-4) '_color_cue.png']);
                
                if SHOW
                    plot_prior_saliency(prior_source, cueImgCnt,cue_by_color, prior_source_with_cue, show_para);
                    figure('name','prior_source_with_cue'); imshow(prior_with_cue_map);
                end
            end
    end
    
    %% second graph
    switch graphOpts.second_graph_method,
        case 'feature'
            aff_mat = gene_weight( spAdjcMat, sp_fea, theta, graphOpts, bd,[],2);
        case 'multi'
            aff_mat = gene_weight( spAdjcMat, [sp_fea  repmat(prior_source_with_cue,1,3)], theta, graphOpts, bd,[],2);
    end
    %% gene foreground seeds %
    pos_label_inds = find(BW_sp==1);
    num_seed = fseedOpts.fseed_num;
    [~,idx] = sort(-prior_source_with_cue);
    im_label_inds = idx(1:num_seed);
    pos_label_inds = intersect(pos_label_inds, im_label_inds);
    
    
    sp_saliency = solve_les_with_dirichlet_v2( aff_mat, pos_label_inds, prior_source_with_cue(pos_label_inds),neg_label_inds, prior_source_with_cue, submodularOpts);
    pos_label_inds_choosen = pos_label_inds;
    
    %% assign the saliency value to each pixel
    [label_im] = draw_seed(input_im, superpixels,  pos_label_inds_choosen, [255,0,0]);
    soutname=[saldir imnames(ii).name(1:end-4) '_seed.png'];
    imwrite(label_im,soutname);
    [im_saliency] = saliency_sp2im(sp_saliency, sp_inds, sp_num, m, n, w);
    
    % post process: just smoothing!
    myfilter = fspecial('gaussian',[5 5], 1);
    im_saliency = imfilter(im_saliency, myfilter, 'replicate');
    
    %% output
    if SHOW
        figure('name', 'saliency map');
        imshow(im_saliency);
         figure('name', 'seed ');
        imshow(label_im);
    end
    outname=[saldir imnames(ii).name(1:end-4) '.bmp'];
    imwrite(im_saliency,outname);
    
end
% test_evaluate