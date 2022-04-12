function [ranks, pos_label_inds_choosen, obj_val] = solve_submodular_naive(...
    aff_mat,seed_num, pos_label_inds_all, neg_label_inds, prior_source, opts,input_im, superpixels,show_para,SHOW,saldir,imname)

%
num_vtx = size(aff_mat, 1);
num_pos_label = length(pos_label_inds_all);
if seed_num > num_pos_label;
    seed_num = num_pos_label;
    disp(['new posotive seed number: ' num2str(seed_num)]);
end

obj_val = [];
current_pos_label_inds = [];
current_ranks = prior_source;
% y_source = prior_source;
seedtype = getoptions(opts,'seedtype','prior');
ytype = getoptions(opts,'ytype','prior');
gamma = getoptions(opts,'gamma','prior');
single_label_saliency = ones(size(prior_source)); %jjcao 


for i= 1:seed_num    
    cand_pos_label_inds = setdiff(pos_label_inds_all, current_pos_label_inds);
    
    num_cand = length(cand_pos_label_inds);
    cand_ranks = zeros(num_vtx, num_cand);
    cand_obj_val = zeros(num_cand, 1);
    for j = 1:num_cand 
        pos_label_inds = [current_pos_label_inds; cand_pos_label_inds(j)];
        
        switch lower(seedtype)
            case 'initial'
                pos_label_ranks = prior_source(pos_label_inds);
            case 'ones'
                pos_label_ranks = ones(length(pos_label_inds), 1);
            case 'current'
                if i>2
                    pos_label_ranks = current_ranks(pos_label_inds);
                else 
                    pos_label_ranks = prior_source(pos_label_inds);   
                end
        end
        switch lower(ytype)
            case 'initial'
                y_source = prior_source;
            case 'current'
                if i>2
                    y_source = current_ranks;
                else 
                    y_source = prior_source;
                end
        end
        
        cand_ranks(:, j) = solve_les_with_dirichlet_v2(...
            aff_mat, pos_label_inds, pos_label_ranks,neg_label_inds, y_source, opts);
        
        % begin jjcao
%         tmp = min(cand_ranks(:, j)); cand_ranks(:, j)=(cand_ranks(:, j)-tmp)/(max(cand_ranks(:, j))-tmp);
        if i == 1
            single_label_saliency(pos_label_inds) = sum(cand_ranks(:, j)); 
        end    
        % end jjcao
        cand_obj_val(j) = gene_obj_val(cand_ranks(:, j), pos_label_inds, gamma,prior_source, single_label_saliency);
    end
    
    [max_val, max_ind] = max(cand_obj_val);
    
    if i > 1 && max_val < obj_val(i - 1)        
%         current_pos_label_inds = [current_pos_label_inds; cand_pos_label_inds(max_ind)];
        break;
    end
    
    current_pos_label_inds = [current_pos_label_inds; cand_pos_label_inds(max_ind)];
    obj_val = [obj_val; max_val];
    current_ranks = cand_ranks(:, max_ind);        
    current_ranks=(current_ranks-min(current_ranks(:)))/(max(current_ranks(:))-min(current_ranks(:)));
     

    if SHOW
        [label_im] = draw_seed(input_im, superpixels, current_pos_label_inds);
        [rank_im] = saliency_sp2im(current_ranks, show_para.sp_inds, show_para.sp_num, show_para.m,show_para.n, show_para.w);

        if i==1
           h = figure; clf 
           h =draw_maps(label_im,rank_im,h);
        else
             h =draw_maps(label_im,rank_im,h);
        end
    end    
    
end

ranks = current_ranks;
pos_label_inds_choosen = current_pos_label_inds;

% if SHOW
% figure;set(gcf,'color','white');
% plot(obj_val); 
% % saveas(gcf,[saldir, imname(1:end-4), '_sub_curve.jpg'])
% end

