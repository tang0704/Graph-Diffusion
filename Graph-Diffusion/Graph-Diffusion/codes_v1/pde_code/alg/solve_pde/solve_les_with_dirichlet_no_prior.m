function [ranks] = solve_les_with_dirichlet_no_prior(aff_mat, label_ranks, label_inds, opts)

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
%   opts : 
%                         solver: 0 -- LU decomposition.  1 -- naive backslash
%                         alpha : parameter for the diffusion sourse
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

if ~exist('opts', 'var')
    opts = [];
end

num_vtx = size(aff_mat, 1); % N
%num_label = length(label_inds);

% convert 'label_ranks' to a column vector
if size(label_ranks, 2) > 1
    label_ranks = label_ranks';
end
if size(label_inds, 2) > 1
    label_inds = label_inds';
end

ground_conductance_coeff = getoptions(opts,'ground_conductance_coeff',0.1);
% consider wheather add a ground conductance vertix to boundary condition for LES 
ground_cond = ground_conductance_coeff.*ones(num_vtx, 1); % tune this parameters
aff_mat = [[aff_mat, ground_cond]; [ones(1,num_vtx), 0]] ;
aff_mat = fn_normout(aff_mat) ;
label_ranks = [label_ranks; 0];
num_vtx_all = num_vtx + 1;
label_inds = [label_inds; num_vtx_all];
ulabel_inds = 1 : num_vtx_all;
ulabel_inds(label_inds) = [];

% Compute Laplace
lap_mat = full(diag(sum(aff_mat,2)) - aff_mat) ;
% lap_mat = lap_mat*lap_mat;
% alpha = getoptions(opts, 'alpha', 0);
% lap_mat = lap_mat + alpha*eye(num_vtx_all);

A = lap_mat(ulabel_inds, ulabel_inds); 
%b = - sum(lap_mat(ulabel_inds, label_inds), 2);
b = - lap_mat(ulabel_inds, label_inds)*label_ranks;

% num_ulabel = length(ulabel_inds);
% alpha = getoptions(opts, 'alpha', 0);
% A = lap_mat(ulabel_inds, ulabel_inds) + alpha*eye(num_ulabel, num_ulabel);
% b = - sum(lap_mat(ulabel_inds, label_inds), 2);

solver = getoptions(opts, 'solver', 0);
ulabel_ranks = solve_axb(A, b, solver);

ranks = zeros(num_vtx, 1);
ranks(label_inds(1:end-1)) = label_ranks(1:end-1); % remove ground conductance 
ranks(ulabel_inds) = ulabel_ranks;

% switch lower(method) 
%     case 'hard'
%         %b = -sum(lap_mat(ulabel_inds, label_inds), 2) ;
%         ulabel_ranks = solve_axb(A, b, solver);
%         ranks = zeros(N, 1);
%         ranks(label_inds) = label_ranks;
%         ranks(ulabel_inds) = ulabel_ranks;
%     case 'soft'
%         for i = length(label_inds)
%             add = zeros(1, N);
%             add(label_inds(i)) = 1;
%             A = [lap_mat; add];
%             b = [b; label_ranks(i)];
%         end
%         ranks = solve_axb(A, b, solver);
%     otherwise
%         error('Unknown method.');
% end

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
