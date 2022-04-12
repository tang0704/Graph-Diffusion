function [ranks, pos_label_inds_choosen, obj_val] = solve_submodular_naive1(...
    aff_mat,s_num, pos_label_inds_all, neg_label_inds, prior_source, opts,input_im, superpixels,show_para,SHOW,saldir,imname,seeds_sp)

%
num_vtx = size(aff_mat, 1);
num_pos_label = length(pos_label_inds_all);


obj_val = [];
current_pos_label_inds = [];
current_ranks = prior_source;
% y_source = prior_source;
seedtype = getoptions(opts,'seedtype','prior');
ytype = getoptions(opts,'ytype','prior');
gamma = getoptions(opts,'gamma','prior');
single_label_saliency = ones(size(prior_source)); %jjcao 



seeds_sp1=[];
% joints=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25];
joints=[1,2,4,6,10,14,16,18,20,21,23,25];
% joints=[1,4,6,10,14,16,18,20,21,23,25];
% joints=[1,4,14,16,18,20,21,23,25];
s_num1=length(joints);
for i=1:25
   if ismember(i,joints)  
       seeds_sp1=[seeds_sp1; seeds_sp(i)];
   end    
end    




% tic;
for i = 1:s_num1 
    cand_pos_label_inds = setdiff(pos_label_inds_all, current_pos_label_inds);
    
    num_cand = length(cand_pos_label_inds);
    cand_ranks = zeros(num_vtx,1);
    cand_obj_val=0;
    
    %for j = 1:num_cand
        pos_label_inds = [current_pos_label_inds; seeds_sp1(i)];
         
        
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
        
        cand_ranks(:, 1) = solve_les_with_dirichlet_v2(...
            aff_mat, pos_label_inds, pos_label_ranks,neg_label_inds, y_source, opts);
        
        if i == 1
            single_label_saliency(pos_label_inds) = sum(cand_ranks(:, 1)); 
        end    
        
        cand_obj_val = gene_obj_val(cand_ranks(:, 1), pos_label_inds, gamma,prior_source, single_label_saliency);
        
        current_pos_label_inds = [current_pos_label_inds; seeds_sp1(i)];
        obj_val = [obj_val; cand_obj_val];
        current_ranks = cand_ranks(:, 1);      
        current_ranks=(current_ranks-min(current_ranks(:)))/(max(current_ranks(:))-min(current_ranks(:)));
    %end    
    
end 

% toc;

ranks = current_ranks;
pos_label_inds_choosen = current_pos_label_inds;

% if SHOW
% figure;set(gcf,'color','white');
% plot(obj_val); 
% % saveas(gcf,[saldir, imname(1:end-4), '_sub_curve.jpg'])
% end

