% Graph diffusion for action recognition task
% modified by Tangli Chu based on
% Demo of Learning-Based Linear Elliptic System
% for Saliency Detection
% by Risheng Liu, Guangyu Zhong, JJCAO
 

clc;clear all;close all;
addpath(genpath('.'));
DEBUG = 1; SHOW = 0;

%% gene dir and file name read videos
database = 'MSRA1000'; % MSRA1000, Berkeley300
database_out = 'MSRA1000';
[imdir, spdir, saldir,gdir, seeddir,rgbdir] = gene_dir(database,database_out);

%% data processing
videoname=dir([imdir '*' 'avi']);%video name
videos_num=length(videoname);
 
% disp('Video data processing...');
% 
% read_video(imdir, videoname);
% 
% disp('Data processing done!');

%% parameters group 3: stable
sp_num_max = 200;    % inital superpixel number (upper bound number) 
% theta = 10; % control the edge weight
theta = 15;
% sigmas = 0.089;
% sigmas=0.06;
sigmas=0.07;
% sigmac = 0.125;
sigmac = 0.02;
% sigmac=0.09;
%% back & front seed paameters, 


priorOpts.priortype = 'convex'; % ground : groundtruth as prior; convex : general convex as prior;  matting: use Pan's mehod to matting the convex (todo, maybe not ok)
priorOpts.enlarge_convex_ratio = 0.1;  % used for enlarging foreground convex
priorOpts.bseed_thresh = 26;     % threshold of background 
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

%%  
% loop for each video 
for jj=1:videos_num
% for jj=1:1y
    v_name=videoname(jj).name;% video name
    
    disp('Video data processing...');
    
    % read video and cut by frame
    read_video(imdir, v_name);
    
    disp('Data processing done00!');
    
    
    frame_dir=[imdir 'img_' v_name(1:end-4) '/'];
    imnames=dir([frame_dir '*' 'bmp']);% read all frames
    
    frame_num=length(imnames);
    
    fid=fopen([seeddir v_name(1:end-8) '.skeleton'], 'r');% read the skeleton file corresponding to the video
    
    whole_frame_num=fgetl(fid);% read first line, which is the number of frame 
    % save rgb info to '.txt' 
    fid1=fopen([rgbdir v_name(1:end-8) '.txt'], 'at+');
    % the fisrt line of '.txt' is the number of frame 
    fprintf(fid1, '%d\n', frame_num);
    
    sprintf('%d: %s', jj, v_name)
    sprintf('%s%d%s%d%s\n','processing the ', jj, '/',videos_num ,'th video...')
    
    % loop for each frame
    for ii=1:frame_num
        close all;    clear w;
        %% read image and preprocessing: : just smoothing!
        imname=imnames(ii).name;% read the name of img
        num_ii=(ii-1)*10+1;% read by order
        imname=[imname(1:25) num2str(num_ii) '.bmp'];
        %sprintf('%d: %s', ii, imname)
        [input_im_frame] = imread([frame_dir imname]);    
        
        [seeds_location1, seeds_location2, per_num,numerror,numless1,data0] = gene_seeds_location(fid,ii,whole_frame_num,v_name);% read 2d coordinates of 25 pixels 
        
        % In return, the number of persons, the number of joints, rgb values into '.txt'
         
        % skip this video if there're all 0 for a row in '.txt'
        if data0==1
            copyfile([imdir v_name], 'D:\DATA\error_video\dirty\');
            copyfile([imdir v_name], 'D:\DATA\rgb008\');
            break;
        end
        
        
        if per_num>0
            fprintf(fid1, '%d\n', per_num);% the number of persons in each frame
            fprintf(fid1, '%d\n', 25);% the number of joints
        end    
        
        
        
        % skip this vedio if there's mistake for the number of frames
        if numerror==1
            copyfile([imdir v_name], 'D:\DATA\error_video\num3\');
            copyfile([imdir v_name], 'D:\DATA\rgb008');
            break;
        end
        
         % kip this vedio if there's mistake for the number of persons
        if numless1==1
            copyfile([imdir v_name], 'D:\DATA\error_video\less1\');
            copyfile([imdir v_name], 'D:\DATA\rgb008');
            break;
        end
        
        
        % crop the img according to the number of persons 
        for p=1:per_num
            % write rgb values for each frame 
            rgb=zeros(25, 3);
            
            if p==1
                seeds_location=seeds_location1;
            else
                seeds_location=seeds_location2;
            end
            
            [input_im,seeds_location] = gene_seeds_location_crop(input_im_frame,seeds_location);
%             crop_imdir=[frame_dir 'crop/'];
            crop_imdir=[imdir 'crop_' v_name(1:end-4) '/'];
            
            if ~isdir(crop_imdir)
                mkdir(crop_imdir);
            end
            imname=[imname(1:end-4) '_crop_p' num2str(p) '.bmp'];
            imwrite(input_im,[crop_imdir imname],'bmp');% save img after cropping
            
            
            [m,n,k] = size(input_im);
            w = [m,n,1,m,1,n];
            
            %     imshow(input_im);
            %     hold on;
            %     x=1;y=1 ;
            %     plot(x,y,'ro');
            %     hold off;
            %%[pixels,sp,in,cen,n]=gene_pixel(imdir,imname,m,n);
            
            %% generate superpixels
            [superpixels, spAdjcMat, sp_inds, sp_center, sp_npix] = gene_superpixel(crop_imdir, imname, sp_num_max, spdir, m, n);
            %     [superpixels, spAdjcMat, sp_inds, sp_center, sp_npix] = run_superpixel(input_im, sp_num_max) ;  
            graphOpts.npix = sp_npix;
            sp_num = size(spAdjcMat,1);% the number of superpixels
            
            
            %% read 25 seeds data (add by T 10-15)
            seeds_sp = gene_sp_seeds(seeds_location, superpixels);
            s_num = size(seeds_sp,1);% the number of seeds, 25
            
            
            %% compute the feature (mean color in lab color space)
            [sp_fea, rgb_fea]  = gene_feature(input_im, superpixels, sp_center, sp_npix, graphOpts);
            min_sp_fea = min(sp_fea(:)); min_rgb_fea = min(rgb_fea(:));
            sp_fea=(sp_fea-min_sp_fea)/(max(sp_fea(:))-min_sp_fea);% lab color values
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
                    % add by T 10-19
                    %[neg_label_inds, EBW_sp, CBW_sp, BW_sp] = gene_convex_seeds(spAdjcMat,seeds_location, superpixels, sp_center,input_im, priorOpts,SHOW,saldir,imname);
                    [neg_label_inds, EBW_sp, CBW_sp, BW_sp] = gene_convex_seeds1(spAdjcMat,seeds_location, superpixels, sp_center,input_im, priorOpts,SHOW,saldir,imname);% Harris+25 
                    
                    
                    cueImgCnt = cue_by_img_center(sp_center, superpixels);%sp_center, sp_npix % only center & color priors be used here
                    tmp = min(cueImgCnt); cueImgCnt=(cueImgCnt-tmp)/(max(cueImgCnt)-tmp);
                    %             [c1, cue_by_color] = cue_by_img_color(train_data_dir, superpixels,sp_center,sp_inds,rgb_fea,show_para);
                    %             tmp = min(cue_by_color); cue_by_color=(cue_by_color-tmp)/(max(cue_by_color)-tmp);
                    
                    bd=unique([superpixels(1,:),superpixels(m,:),superpixels(:,1)',superpixels(:,n)']); 
                    bd  = setdiff(bd',neg_label_inds);
                    bd = [neg_label_inds;bd]; % for smoothing 
                    
%                     aff_mat = gene_weight( spAdjcMat, sp_fea, theta, graphOpts,bd,[],2);% affinity matrix,
                    aff_mat = gene_weight( spAdjcMat, sp_fea, sigmac, sigmas, graphOpts,bd,[],2,sp_center,[],[]);% affinity matrix,
                    
                    [prior_source] = gene_prior_source(aff_mat, neg_label_inds, priorOpts.doBinary, graphOpts);
                    tmp = min(prior_source); prior_source=(prior_source-tmp)/(max(prior_source)-tmp);
                    
                    
                    %with center prior(add by tang 10-15)
                    prior_source_with_cue =prior_source.*cueImgCnt;% rank
% prior_source_with_cue =prior_source;
                    minpc = min(prior_source_with_cue(:));
                    prior_source_with_cue=(prior_source_with_cue-minpc)/(max(prior_source_with_cue)-minpc);
                    prior_with_cue_map = saliency_sp2im(prior_source_with_cue, ...
                        show_para.sp_inds, show_para.sp_num, show_para.m, ...
                        show_para.n, show_para.w);
%                     imwrite(prior_with_cue_map,[saldir imnames(ii).name(1:end-4) '_p_cue.png']);% generate img 
%                     SHOW=1;
                    if SHOW
                        prior_map = saliency_sp2im(prior_source, ...
                            show_para.sp_inds, show_para.sp_num, show_para.m, ...
                            show_para.n, show_para.w);
                        imwrite(prior_map,[saldir imnames(ii).name(1:end-4) '_background_cue.png']);
                        imshow(prior_map);
                        
                        cueImgCnt_map = saliency_sp2im(cueImgCnt, ...
                            show_para.sp_inds, show_para.sp_num, show_para.m, ...
                            show_para.n, show_para.w);
                        imwrite(cueImgCnt_map,[saldir imnames(ii).name(1:end-4) '_center_cue.png']);
                    end
                    
                    
                    %annotated by T 10-10
                    %             prior_source_with_cue =prior_source.*cueImgCnt.*cue_by_color;
                    %             minpc = min(prior_source_with_cue(:));
                    %             prior_source_with_cue=(prior_source_with_cue-minpc)/(max(prior_source_with_cue)-minpc);
                    %
                    %             prior_with_cue_map = saliency_sp2im(prior_source_with_cue, ...
                    %                 show_para.sp_inds, show_para.sp_num, show_para.m, ...
                    %                 show_para.n, show_para.w);
                    %             imwrite(prior_with_cue_map,[saldir imnames(ii).name(1:end-4) '_p_cue.png']);%
                    %
                    %
                    %             if SHOW
                    %                 prior_map = saliency_sp2im(prior_source, ...
                    %                     show_para.sp_inds, show_para.sp_num, show_para.m, ...
                    %                     show_para.n, show_para.w);
                    %                 imwrite(prior_map,[saldir imnames(ii).name(1:end-4) '_background_cue.png']);
                    %
                    %                 cueImgCnt_map = saliency_sp2im(cueImgCnt, ...
                    %                     show_para.sp_inds, show_para.sp_num, show_para.m, ...
                    %                     show_para.n, show_para.w);
                    %                 imwrite(cueImgCnt_map,[saldir imnames(ii).name(1:end-4) '_center_cue.png']);
                    %
                    %                 cue_by_color_map = saliency_sp2im(cue_by_color, ...
                    %                     show_para.sp_inds, show_para.sp_num, show_para.m, ...
                    %                     show_para.n, show_para.w);
                    %                 imwrite(cue_by_color_map,[saldir imnames(ii).name(1:end-4) '_color_cue.png']);
                    %
                    %                 if SHOW
                    %                     plot_prior_saliency(prior_source, cueImgCnt,cue_by_color, prior_source_with_cue, show_para);
                    %                     figure('name','prior_source_with_cue'); imshow(prior_with_cue_map);
                    %                 end
                    %             end
            end
            
            %% gene foreground seeds
            pos_label_inds = find(BW_sp);
            
            
            %% perform learning-based les / select seeds!  
            switch graphOpts.second_graph_method,
                case 'feature'
%                     aff_mat = gene_weight( spAdjcMat, sp_fea, theta, graphOpts, bd,[],2);
                    aff_mat = gene_weight( spAdjcMat, sp_fea, sigmac, sigmas, graphOpts, bd, seeds_sp,2, sp_center,neg_label_inds,pos_label_inds);
%                       aff_mat = gene_weight( spAdjcMat, sp_fea, sigmac, sigmas, graphOpts, bd, seeds_sp,3, sp_center,[],[],[]);
                case 'multi'
                    aff_mat = gene_weight( spAdjcMat, [sp_fea  repmat(prior_source_with_cue,1,3)], theta, graphOpts, bd,[],2);
            end
            
            seed_num = numel(pos_label_inds);% the number of seeds which belong to foreground 
            
            %annotated by T 10-10
%                 [sp_saliency, pos_label_inds_choosen] = solve_submodular_naive(...
%                     aff_mat, seed_num, pos_label_inds, neg_label_inds, prior_source_with_cue, ...
%                     submodularOpts,input_im,...
%                     superpixels,show_para,SHOW,saldir,imname); % select heat sources automaticly
%             
            [sp_saliency, pos_label_inds_choosen] = solve_submodular_naive1(...
                aff_mat, s_num, pos_label_inds, neg_label_inds, prior_source_with_cue, ...
                submodularOpts,input_im, superpixels,show_para,SHOW,saldir,imname, seeds_sp); % consider 25 seeds as heat sources 
            
            
            %% assign the saliency value to each pixel
            [label_im] = draw_seed(input_im, superpixels,  pos_label_inds_choosen, [255,0,0]);
            soutname=[saldir imnames(ii).name(1:end-4) '_seed.png'];
            imwrite(label_im,soutname);
            [im_saliency] = saliency_sp2im(sp_saliency, sp_inds, sp_num, m, n, w);
            
            
%             [label_im1] = draw_seed(input_im, superpixels,  pos_label_inds, [255,0,0]);
%             soutname1=[saldir imnames(ii).name(1:end-4) '_ffff.png'];
%             imwrite(label_im1,soutname1);
%             
%             [label_im2] = draw_seed(input_im, superpixels,  neg_label_inds, [255,0,0]);
%             soutname2=[saldir imnames(ii).name(1:end-4) '_bbbb.png'];
%             imwrite(label_im2,soutname2);
            
            
            
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
            % save saliency map
            outname=[saldir imnames(ii).name(1:end-4) '.bmp'];
            imwrite(im_saliency,outname);
            
            
            rgb=gene_RGB(s_num, spAdjcMat, seeds_sp, superpixels, input_im, sp_saliency);
            
            for seed=1:s_num
               fprintf(fid1,'%d ',floor(rgb(seed,1))); 
               fprintf(fid1,'%d ',floor(rgb(seed,2))); 
               fprintf(fid1,'%d\n',floor(rgb(seed,3))); 
            end
            
            
            if per_num==2 && p==1
                fprintf(fid1,'%d\n', 25); 
            end    
            
        end
    end
    
    fclose(fid);
    fclose(fid1);
    rmdir(frame_dir,'s');
    if isdir(spdir)
        rmdir(spdir,'s');
    end
    
    if numerror==1 || numless1==1 || data0==1
       delete([rgbdir v_name(1:end-8) '.txt']);
       continue;
    end    
    % skip if there's 3 persons in one frame
    rmdir(crop_imdir, 's');
end
sprintf('%s', 'LTD done!') 
% test_evaluate
