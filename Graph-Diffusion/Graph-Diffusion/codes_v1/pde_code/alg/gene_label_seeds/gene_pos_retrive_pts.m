function  label_inds = gene_pos_retrive_pts(prior,CBW_sp, opts)
% sort the prior    the most num salient seeds are foreground seeds
% GuangyuZhong, JJCAO, 11/10/2013

opts.null = [];
fseed_percentage = getoptions(opts,'fseed_percentage', 1);
fseed_sort = getoptions(opts,'fseed_sort', true);

[inds, ~] = find(CBW_sp);
if fseed_percentage == 1;
    label_inds = inds;
else
    num = ceil(fseed_percentage*numel(inds));
    if fseed_sort
        [~,idx] = sort(-prior);
        im_label_inds = idx(1:num);
        label_inds = intersect(inds, im_label_inds);    
    else
        idx = randperm(numel(inds));
        label_inds = inds(idx(1:num));
    end
end

