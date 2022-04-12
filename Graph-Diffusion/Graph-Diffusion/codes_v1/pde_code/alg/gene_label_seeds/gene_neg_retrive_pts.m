function  [label_inds, EBW_sp, CBW_sp, BW_sp] = gene_neg_retrive_pts(spAdjcMat, superpixels, input_im, opt,DEBUG,SHOW,saldir,imname)
% use convex method to get background seeds
% Guangyuzhong 11/10/2013
%
% jjcao @ 2013
%

%% compute convex
bseed_thresh = getoptions(opt,'bseed_thresh',26);
enlarge_convex_ratio = getoptions(opt,'enlarge_convex_ratio',0.2);
priortype = getoptions(opt,'priortype','convex');

[corner_im2,EnIm] = getsalientpoints(input_im, enlarge_convex_ratio);
% if strcmp(priortype, 'matting')
%     corner_im = corner_im2;
% else
% corner_im = elimatepoint(corner_im2,bseed_thresh, superpixels, spAdjcMat);
corner_im = elimatepoint(corner_im2,bseed_thresh, [], []);
% end

[row,col] = size(corner_im);
[y,x] = ind2sub([row,col],find(corner_im == 1));%find positiion of the corner points
dt = DelaunayTri(x,y);
if(~size(dt,1))
    return;
end

%% obtain the convex hull
[K, AV] = convexHull(dt);%find the points to plot the convex hull
% if AV/(row*col) > 0.6
%     corner_im = elimatepoint(corner_im2,bseed_thresh, [], []);
%     [y,x] = ind2sub([row,col],find(corner_im == 1));%find positiion of the corner points
%     dt = DelaunayTri(x,y);
%     if(~size(dt,1))
%         return;
%     end
%     [K, AV] = convexHull(dt);
% end

BW = roipoly(corner_im,x(K),y(K));%obtain the pixels inside the convex hull

%% matting, if needed. Actually it is bad.
if strcmp(priortype, 'matting')
    if max(double(input_im(:)))>12;
        input_im = double(input_im)/255;
    end
    r_guide = 120;
    eps_guide = 10^-6;
    BW_tmp = guidedfilter_color(input_im, double(BW), r_guide, eps_guide);
    level_bw=graythresh(BW_tmp);
    BW_new=im2bw(BW_tmp,level_bw);
    BW_new = elimatepoint(BW_new,bseed_thresh, [], []);
    BW_new = imfill(BW_new,'hole');
    sum(BW(:)) - sum(BW_new(:))
    
    BW = BW_new;
end

%% mark superpixel, in convex hull = 1(foreground), else = 0 (backgound)
sp_num = max(superpixels(:));
BW_sp = zeros(sp_num,1);
for ii = 1:sp_num
    indsp = find(superpixels==ii);
    indsp1 = find(BW(indsp) == 1);
    if (numel(indsp1)>0.2*numel(indsp))
        BW_sp(ii)= 1;
    end
end

%% extend & contract the convex hull
EBW_sp = BW_sp;
tmp = find(BW_sp == 1);
for i = tmp'
    EBW_sp(spAdjcMat(i,:)>0) = 1;
end
if length(find(EBW_sp==1))/sp_num > 0.9
    EBW_sp = BW_sp;
end

CBW_sp = BW_sp;
if ~strcmp(priortype, 'matting')
    tmp = find(BW_sp == 0);
    for i = tmp'
        CBW_sp(spAdjcMat(i,:)>0) = 0;
    end
end

if SHOW        
    figure; subplot(2,2,1);imshow(EnIm);hold on;
    [y2,x2] = ind2sub(size(input_im),find(corner_im2 ~= 0));
    scatter(x2,y2,50,[1 0 0],'filled'); 
    plot(dt.X(K,1),dt.X(K,2), 'r');
    
    subplot(2,2,2);
    [Elabel_im_neg] = draw_seed(input_im, superpixels, find(EBW_sp==1));    
    imshow(Elabel_im_neg); hold on;    
    scatter(x,y,50,[0, 0 1],'filled'); 
    plot(dt.X(K,1),dt.X(K,2), 'g');
   
    subplot(2,2,3);
    [label_im_neg] = draw_seed(input_im, superpixels, find(BW_sp==1));    
    imshow(label_im_neg); hold on;    
    plot(dt.X(K,1),dt.X(K,2), 'g');

    subplot(2,2,4);
    Clabel_im_neg = draw_seed(input_im, superpixels, find(CBW_sp==1));    
    imshow(Clabel_im_neg); hold on;    
    plot(dt.X(K,1),dt.X(K,2), 'g');

end   

if DEBUG        
    figure('name','Enimage');set(gcf,'color','white');imshow(EnIm); 
    EnImname = [saldir, imname(1:end-4), '_En','.jpg'];
    imwrite(EnIm,EnImname);hold on;
    
    [y2,x2] = ind2sub(size(input_im),find(corner_im2 ~= 0));
    scatter(x2,y2,50,[1 0 0],'filled'); 
    plot(dt.X(K,1),dt.X(K,2), 'r');
    harrisname = [saldir, imname(1:end-4), '_harris','.jpg'];
    saveas(gcf, harrisname);
    
    figure('name','enlarge');set(gcf,'color','white');
    [Elabel_im_neg] = draw_seed(input_im, superpixels, find(EBW_sp==1));    
    imshow(Elabel_im_neg); hold on;
    Elabel_im_neg_name = [saldir, imname(1:end-4), '_enlarge','.jpg'];
    imwrite(Elabel_im_neg,Elabel_im_neg_name);
    
    scatter(x,y,50,[0, 0 1],'filled'); 
    plot(dt.X(K,1),dt.X(K,2), 'g');
    saveas(gcf,[saldir, imname(1:end-4), '_enlarge_point','.jpg']);
    
    figure('name','initial_neg');set(gcf,'color','white');
    [label_im_neg] = draw_seed(input_im, superpixels, find(BW_sp==1));    
    imshow(label_im_neg); hold on; 
    plot(dt.X(K,1),dt.X(K,2), 'g');
    saveas(gcf,[saldir, imname(1:end-4), '_label_neg','.jpg']);            

    figure('name','Contrac');set(gcf,'color','white');
    Clabel_im_neg = draw_seed(input_im, superpixels, find(CBW_sp==1));  
    imshow(Clabel_im_neg); hold on;    
    plot(dt.X(K,1),dt.X(K,2), 'g');
    saveas(gcf,[saldir, imname(1:end-4), '_clabel_neg','.jpg']);

end

%% select background/negative retrive points
bseed_percentage = getoptions(opt,'bseed_percentage',1);
[label_inds] = find(EBW_sp==0 ); % find points outside the convex

if bseed_percentage < 1
    num_all = length(label_inds);
    num = floor(bseed_percentage*num_all);    

    [EnIm_sp] = img2superpixel(EnIm,superpixels);
    [~,idx] = sort(EnIm_sp(label_inds));
    idx = idx(1:num);
    
    label_inds = label_inds(idx);
end

