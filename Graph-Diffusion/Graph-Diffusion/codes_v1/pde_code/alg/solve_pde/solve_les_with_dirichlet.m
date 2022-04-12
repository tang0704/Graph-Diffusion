function [ranks] = solve_les_with_dirichlet(aff_mat, pos_label_ranks, pos_label_inds, ...
    neg_label_inds, prior_source, opts)

%
%	Solve linear elliptic system (LES) with dirichlet boundary for semi-supervised
%	ranking
%
% -------
% INPUTS:
%
%	aff_mat        : affinity matrix ([N x N]) 
%                     where N is the number of vertices
%   label_ranks     : rankings for labeled points
%   label_inds       : indexs for labeled points
%   prior_source   : prior source function of unlabled points for LES
%   opts : 
%                         solver: 0 -- LU decomposition.  1 -- naive backslash
%                         alpha : parameter for the diffusion sourse
%                         beta: parameter for the 'prior_source'
%                         add_ground: 0 -- do not add ground conductance vertix,
%                         1 -- add ground conductance vertix.
%
% --------
% OUTPUTS:           
%   ranks       : rankings for  points
%
% -------------
%
% Written by Risheng Liu @ DUT, 2013.
% All rights reserved.
%

% label_ranks = unique(label_ranks);
% num_label = length(label_ranks);

% Compute Laplace
% lap_mat = full(diag(sum(aff_mat,2)) - aff_mat) ;
% lap_mat = aff_mat;


if ~exist('opts', 'var')
    opts = [];
end

num_vtx = size(aff_mat, 1); % N
%num_label = length(label_inds);

% convert 'label_ranks' to a column vector
if size(pos_label_ranks, 1) == 1
    pos_label_ranks = pos_label_ranks';
end
if size(pos_label_inds, 1) == 1
    pos_label_inds = pos_label_inds';
end

% consider wheather add a ground conductance vertix to boundary condition for LES 
ground_cond = 0.01*ones(num_vtx, 1); % tune this parameters
ground_cond(neg_label_inds) = 1;
ground_cond(pos_label_inds) = 0;

aff_mat = [[aff_mat, ground_cond]; [ones(1,num_vtx), 0]] ;
%aff_mat = [[aff_mat, ground_cond]; ground_cond', 0];
aff_mat = fn_normout(aff_mat);
label_ranks = [pos_label_ranks; 0];
num_vtx_all = num_vtx + 1;
label_inds = [pos_label_inds; num_vtx_all];
ulabel_inds = 1 : num_vtx_all;
ulabel_inds(label_inds) = [];

alpha = getoptions(opts, 'alpha', 0);
% Compute Laplace
lap_mat = full(diag(sum(aff_mat,2)) - aff_mat);
lap_mat = lap_mat + alpha*eye(num_vtx_all);

A = lap_mat(ulabel_inds, ulabel_inds);
%b = - sum(lap_mat(ulabel_inds, label_inds), 2);
b = - lap_mat(ulabel_inds, label_inds)*label_ranks;
label_ranks

b = b + alpha*prior_source(ulabel_inds);

solver = getoptions(opts, 'solver', 0);
ulabel_ranks = solve_axb(A, b, solver);

%ulabel_ranks=(ulabel_ranks-min(ulabel_ranks(:)))/(max(ulabel_ranks(:))-min(ulabel_ranks(:)));

ranks = zeros(num_vtx, 1);
ranks(label_inds(1:end-1)) = label_ranks(1:end-1); % remove ground conductance 
ranks(ulabel_inds) = ulabel_ranks;

%ranks=(ranks-min(ranks(:)))/(max(ranks(:))-min(ranks(:)));

function [x] = solve_axb(A, b, solver)

switch solver 
    case 0,
        if det(A) == 0
            ilu_setup = [];
            ilu_setup.type = 'ilutp';
            [L, U, P] = ilu(sparse(A), ilu_setup);
            [p, ~] = find(P);  p(p) = 1:length(p) ;
            x = U\(L\(b(p,:)));
        else
            x = A \ b;
        end
    case 1
        x = A \ b;
    otherwise
        error('Unknown solver.');
end
