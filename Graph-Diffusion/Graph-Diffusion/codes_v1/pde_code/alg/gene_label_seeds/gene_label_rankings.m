function [pos_label_inds, neg_label_inds] = gene_label_rankings(superpixels, input_im, pos_num, neg_num, opt)
%pos_label_inds    points in the convex hull
%pos_label_vals     = 1      [pos_num x 1]
%neg_label_inds    points in background
%neg_label_vals     = 0      [neg_num x 1]
% fb  = 1  choose foreground seeds  
% fb  = 0  choose background seeds
% Guangyuzhong 26/09/2013
% modified by Risheng Liu
% divide pos and neg functions  by  Guangyuzhong 09/10/2013
% too many bugs!!!
thresh_edge = getoptions(opt,'thresh_edge',50);
thresh = getoptions(opt,'thresh',26);
pos_type = getoptions(opt,'pos_type','convex');
neg_type = getoptions(opt,'neg_type','convex');
pind = getoptions(opt,'pind','sort');
nind = getoptions(opt,'nind','sort');

[pos_label_inds] = gene_rankings(superpixels, input_im, pos_num, pos_type, 1,pind,thresh,thresh_edge);
[neg_label_inds] = gene_rankings(superpixels, input_im, neg_num, neg_type, 0,nind,thresh,thresh_edge);



function  [label_inds] = gene_rankings(superpixels, input_im, num, type, fb,ind,thresh,thresh_edge)
        if ~exist('thresh','var')
            thresh = 26;
        end
        [corner_im2,EnIm] = getsalientpoints(input_im);        
        corner_im = elimatepoint(corner_im2,thresh);
        [row,col] = size(corner_im);
        [y,x] = ind2sub([row,col],find(corner_im == 1));%find positiion of the corner points
        dt = DelaunayTri(x,y);
        if(~size(dt,1))
            return;
        end
        [EnIm_sp] = img2superpixel(EnIm,superpixels); 
        % obtain the convex hull
        [k, av] = convexHull(dt);%find the points to plot the convex hull
        BW = roipoly(corner_im,x(k),y(k));%obtain the pixels inside the convex hull

switch lower(type)

        
    case 'matting'
         if ~exist('thresh_edge','var')
            thresh_edge = 50;
        end
        if max(double(input_im(:)))>12;
            input_im = double(input_im)/255;
        end
        r_guide = 60;
        eps_guide = 10^-6;
        BW_tmp = guidedfilter_color(input_im, double(BW), r_guide, eps_guide);
        level_bw=graythresh(BW_tmp);
        BW_new=im2bw(BW_tmp,level_bw);
        BW_new = elimatepoint(BW_new,thresh_edge);
        BW_new = imfill(BW_new,'hole');
        BW = BW_new;
        
end
 [BW_sp] = img2superpixel(BW,superpixels);

 [label_inds, label_vals] = find(BW_sp==fb);
  num_all = length(label_inds);
if num > num_all
    num = num_all;
    if fb
    disp(['new posotive seed number: ' num2str(num)]);
    else
      disp(['new negtive seed number: ' num2str(num)]);  
    end
end

switch lower(ind)
    case 'random'
        idx = randperm(num_all);
        idx = idx(1:num);
    case 'sort'
       if fb
        [~,idx] = sort(-EnIm_sp(label_inds));
        idx = idx(1:num);
       else
        [~,idx] = sort(EnIm_sp(label_inds));
        idx = idx(1:num);  
       end
end

label_inds = label_inds(idx);
label_vals = label_vals(idx);


