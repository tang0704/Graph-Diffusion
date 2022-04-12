function plot_prior_saliency(ranks, cueImgCnt, cue_by_color, ranks_with_prior, show_para)

figure('name','original, center and merged prior');

subplot(2,2,1);
[im_saliency] = saliency_sp2im(ranks, ...
    show_para.sp_inds, show_para.sp_num, show_para.m, show_para.n, show_para.w);
imshow(im_saliency);


subplot(2,2,2);
[im_saliency] = saliency_sp2im(cueImgCnt, ...
    show_para.sp_inds, show_para.sp_num, show_para.m, show_para.n, show_para.w);
imshow(im_saliency);


subplot(2,2,3);
[im_saliency] = saliency_sp2im(cue_by_color, ...
    show_para.sp_inds, show_para.sp_num, show_para.m, show_para.n, show_para.w);
imshow(im_saliency);

subplot(2,2,4);
[im_saliency] = saliency_sp2im(ranks_with_prior, ...
    show_para.sp_inds, show_para.sp_num, show_para.m, show_para.n, show_para.w);
imshow(im_saliency);
