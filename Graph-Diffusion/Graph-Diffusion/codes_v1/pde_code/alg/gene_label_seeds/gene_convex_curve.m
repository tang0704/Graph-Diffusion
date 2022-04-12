function  [label_inds,EBW_sp, CBW_sp, BW_sp] = gene_convex(spAdjcMat, superpixels, sp_center,input_im, opt,SHOW,saldir,imname)
% use convex method to draw convex (background and foreground)
% Guangyuzhong 11/10/2013
%
% GuangyuZhong @ 04/11/2013
%

%% compute convex
bseed_thresh = getoptions(opt,'bseed_thresh',26);
enlarge_convex_ratio = getoptions(opt,'enlarge_convex_ratio',0.2);
priortype = getoptions(opt,'priortype','convex');

[corner_im2,EnIm] = getsalientpoints(input_im, enlarge_convex_ratio);%T£º¼ÆËãHarris½Çµã
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
    if (numel(indsp1)>0.5*numel(indsp))
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
    [img] = draw_seed(input_im, superpixels, []);
  figure;
    imshow(img); hold on;
%     scatter(x(K),y(K),50,[0, 0 0],'filled');
    plot(dt.X(K,1),dt.X(K,2),'r' ,'LineWidth',2.5);
    saveas(gcf,[saldir, imname(1:end-4), '_foreground_convex.png']);
%     close ;
    figure; imshow(img);hold on;
    pos_label_inds = find(EBW_sp);
    p = sp_center(pos_label_inds,:);
    dtt = DelaunayTri(p(:,2), p(:,1));
    [KK, ~] = convexHull(dtt);
    plot(dtt.X(KK,1),dtt.X(KK,2),'y','LineWidth',2.5);
    saveas(gcf,[saldir, imname(1:end-4), '_background_convex.png']);
%     close;
    
%       pos_label_inds = find(CBW_sp);
%     p = sp_center(pos_label_inds,:);
%     dtt = DelaunayTri(p(:,2), p(:,1));
%     [KK, ~] = convexHull(dtt);
%     hold on; plot(dtt.X(KK,1),dtt.X(KK,2),'g','LineWidth',2.5)
    
%     figure;
%     [Elabel_im_neg] = draw_seed(input_im, superpixels, find(EBW_sp==1));
%     imshow(Elabel_im_neg); hold on;
%     scatter(x,y,50,[0, 0 1],'filled');
    figure; imshow(img); hold on;
    plot(dt.X(K,1),dt.X(K,2),'r' ,'LineWidth',2.5);
    plot(dtt.X(KK,1),dtt.X(KK,2),'y','LineWidth',2.5);
    saveas(gcf,[saldir, imname(1:end-4), '_fore_and_back_convex.png']);
%     close;
end
[EnIm_sp] = img2superpixel(EnIm,superpixels);
[label_inds] = find(EBW_sp==0 ); % find points outside the convex 
% [label_inds] = find(BW_sp==0 & EnIm_sp<min(2*level_mask, 0.9)); 
% [label_inds] = find(BW_sp==0 & EnIm_sp<level_mask); 
bseed_sort = getoptions(opt,'bseed_sort', true);
bseed_percentage = getoptions(opt,'bseed_percentage',0.8);
num_all = length(label_inds);
num = floor(bseed_percentage*num_all);

if bseed_sort
    [~,idx] = sort(EnIm_sp(label_inds));
    idx = idx(1:num);
else
    idx = randperm(num_all);
    idx = idx(1:num);
end

label_inds = label_inds(idx);