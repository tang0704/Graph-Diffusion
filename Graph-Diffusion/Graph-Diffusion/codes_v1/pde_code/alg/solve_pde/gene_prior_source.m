function [prior_source] = gene_prior_source(aff_mat, label_inds, doBinary, opts)

num_label = length(label_inds);
label_ranks = ones(num_label, 1);

prior_source = solve_les_with_dirichlet_no_prior(aff_mat, label_ranks, label_inds, opts);
prior_source=(prior_source-min(prior_source(:)))/(max(prior_source(:))-min(prior_source(:)));
prior_source = 1 - prior_source;

if doBinary
    % generate binary prior map
    th=mean(prior_source);
    prior_source(prior_source<th)=0;
    prior_source(prior_source>=th)=1;
%     prior_source = prior_source_bw'.*prior_source;
end
