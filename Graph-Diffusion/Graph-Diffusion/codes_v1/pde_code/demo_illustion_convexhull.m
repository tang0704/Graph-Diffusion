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
%% back & front seed parameters, 用缩进表示稳定的参数
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
    
    [neg_label_inds, EBW_sp, CBW_sp, BW_sp] = gene_convex_curve(spAdjcMat, superpixels, sp_center,input_im, priorOpts,SHOW,saldir,imname);  
end
% test_evaluate