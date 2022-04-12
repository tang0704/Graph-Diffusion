function [label_inds, EBW_sp, CBW_sp, BW_sp] = gene_convex_seeds1(spAdjcMat,seeds_location ,superpixels, sp_center,input_im, opt,SHOW,saldir,imname)
% use 25 seeds and Harris to draw convex hull C (foreground)
%extend C as C' (outside C', background)

%% convex is from seeds_location
bseed_thresh = getoptions(opt,'bseed_thresh',26);
enlarge_convex_ratio = getoptions(opt,'enlarge_convex_ratio',0.2);
priortype = getoptions(opt,'priortype','convex');

[corner_im2,EnIm] = getsalientpoints(input_im, enlarge_convex_ratio);


headfoot_points = add_headfoot_hull(seeds_location,input_im);

corner_im = elimatepoint(corner_im2,bseed_thresh, [], []);
[row,col] = size(corner_im);
[yy,xx] = ind2sub([row,col],find(corner_im == 1));%find positiion of the corner points

% imshow(input_im);
%     hold on;
%     for i=1:length(xx)
%         x1=xx(i);
%         y1=yy(i);
%         plot(x1,y1,'ro','MarkerFaceColor','r');
%     end    
%     
%     hold off;
% saveas(gcf,[saldir, imname(1:end-4), 'harris.png']);



headfoot_points=[headfoot_points; xx,yy];
% headfoot_points=[xx,yy];

[x]=headfoot_points(:,1);
[y]=headfoot_points(:,2);

% [x]=seeds_location(:,1);
% [y]=seeds_location(:,2);



imshow(input_im);
    hold on;
    for i=1:length(xx)
        x1=x(i);
        y1=y(i);
        plot(x1,y1,'ro','MarkerFaceColor','r');
    end    
    
    hold off;
saveas(gcf,[saldir, imname(1:end-4), 'harris2.png']);




% imshow(input_im);
%     hold on;
%     for i=1:length(x)
%         x1=headfoot_points(i,1);
%         y1=headfoot_points(i,2);
%         plot(x1,y1,'ro');
%     end    
%     
%     hold off;
% saveas(gcf,[saldir, imname(1:end-4), 'all.png']);





dt = DelaunayTri(x,y);

if(~size(dt,1))
    return;
end

%% obtain the convex hull
[K, AV] = convexHull(dt);%find the points to plot the convex hull

BW = roipoly(input_im,x(K),y(K));


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

BBW_sp=EBW_sp;
pos_label_BW_inds=find(BBW_sp);
pp=sp_center(pos_label_BW_inds,:);
ddt=DelaunayTri(pp(:,2), pp(:,1));
[K1, ~] = convexHull(ddt);
% 
% tmp = find(BW_sp == 1);
% for i = tmp'
%     EBW_sp(spAdjcMat(i,:)>0) = 1;
% end
% if length(find(EBW_sp==1))/sp_num > 0.9
%     EBW_sp = BW_sp;
% end




CBW_sp = BW_sp;
if ~strcmp(priortype, 'matting')
    tmp = find(BW_sp == 0);
    for i = tmp'
        CBW_sp(spAdjcMat(i,:)>0) = 0;
    end
end


 SHOW=1;
if SHOW        
    [img] = draw_seed(input_im, superpixels, []);
  figure;
    imshow(img); hold on;

%     plot(ddt.X(K1,1),ddt.X(K1,2),'r' ,'LineWidth',2.5);
%     saveas(gcf,[saldir, imname(1:end-4), '_foreground_convex.png']);
% %     close ;
%     figure; imshow(img);hold on;
%     pos_label_inds = find(EBW_sp);
%     p = sp_center(pos_label_inds,:);
%     dtt = DelaunayTri(p(:,2), p(:,1));
%     [KK, ~] = convexHull(dtt);
%     plot(dtt.X(KK,1),dtt.X(KK,2),'y','LineWidth',2.5);
%     saveas(gcf,[saldir, imname(1:end-4), '_background_convex.png']);

plot(dt.X(K,1),dt.X(K,2),'r' ,'LineWidth',2.5);
    saveas(gcf,[saldir, imname(1:end-4), '_foreground_convex.png']);
%     close ;
    figure; imshow(img);hold on;
%     pos_label_inds = find(EBW_sp);
%     p = sp_center(pos_label_inds,:);
    [K1, ~] = convexHull(ddt);
    plot(ddt.X(K1,1),ddt.X(K1,2),'y','LineWidth',2.5);
    saveas(gcf,[saldir, imname(1:end-4), '_background_convex.png']);

    figure; imshow(img); hold on;
    plot(dt.X(K,1),dt.X(K,2),'r' ,'LineWidth',2.5);
    plot(ddt.X(K1,1),ddt.X(K1,2),'y','LineWidth',2.5);
    saveas(gcf,[saldir, imname(1:end-4), '_fore_and_back_convex.png']);
%     close;
end

[EnIm_sp] = img2superpixel(EnIm,superpixels);
[label_inds] = find(EBW_sp==0 ); % find points outside the convex 
bseed_sort = getoptions(opt,'bseed_sort', false);% modified by T from true to false 10-19
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






