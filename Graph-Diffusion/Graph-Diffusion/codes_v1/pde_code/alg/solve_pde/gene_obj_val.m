function obj_val = gene_obj_val(ranks, label_inds, gamma,prior_source, single_label_saliency)

if ~exist('gamma', 'var')
    gamma = 0.5;
end
if strcmp(gamma,'prior')
%     obj_val = sum(ranks) - sum(1./prior_source(label_inds).^2);
    obj_val = sum(ranks) - sum(1./(1+prior_source(label_inds).^2));
%     obj_val = sum(ranks) - 40*sum(1./prior_source(label_inds).^1);    
else
    %obj_val = sum(ranks(ulabel_inds)) - gamma*sum(ranks(label_inds));
%     obj_val = sum(ranks)- gamma*sum(ones(length(label_inds),1));

%     tmp = sum(1./prior_source(label_inds).^2);
%     tmp = sum(1-prior_source(label_inds));
%     tmp = 1./(single_label_saliency(label_inds));
    tmp = 1./(single_label_saliency(label_inds).*prior_source(label_inds).^2);
%     tmp = length(label_inds);
    tmp = gamma*tmp;
%     tmp(tmp<1) = 1;
    
    obj_val = sum(ranks)- sum(tmp);
end