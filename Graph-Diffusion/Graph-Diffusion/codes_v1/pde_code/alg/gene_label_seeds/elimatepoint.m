function corner_im2 = elimatepoint(corner_im2,thresh, superpixels, aff_mat)
%% function corner_im2 = elimatepoint(corner_im2,thresh)
%%compute salient points
%%Input:
%       corner_im2    : deteced salient points 
%%Output:
%       corner_im2    : salient points after elimating points too close to
%                       the image boundary
%       thresh        : the threshold defining the points too close to the
%                       image boundary

[row,col] = size(corner_im2);
ind = find(corner_im2 == 1);
[x,y] = ind2sub([row,col],ind);
if isempty(aff_mat)
    elipos = [];    
    elipos_x  = find( ( x>(row - thresh) ) | ( x < (thresh + 1))); 
    elipos = [elipos (y(elipos_x)-1)*row + x(elipos_x)];

    elipos_y  = find( (y > (col - thresh))  |  (y < (thresh + 1)));
    elipos = [elipos; (y(elipos_y)-1)*row + x(elipos_y)];

    corner_im2 (elipos)=0; 
else
    n_ind = filter_by_sp_neighbors(ind, superpixels, aff_mat);
    corner_im2(n_ind)=0; 
end

function n_ind = filter_by_sp_neighbors(ind, superpixels, aff_mat)
n_ind = ind;
sp = superpixels(ind);
for i = 1:length(sp)
%     neighsp = sp(i);
    neighsp = find(aff_mat(sp(i),:)~=0);
    
    tmpSp = sp;
    tmpSp(i) = [];
    C = intersect(tmpSp, neighsp);
    
    if ~isempty(C)
        n_ind(i) = -1;        
    end
end
n_ind(n_ind == -1) = [];

% function elipos_x = filter_by_sp_neighbors(elipos_x, x, y, superpixels, aff_mat)
% sp_x = superpixels(x(elipos_x), y(elipos_x));
% for i = 1:length(sp_x)
%     neighsp = find(aff_mat(sp_x(i),:)~=0);
%     tmpx = x;
%     tmpx(elipos_x(i)) = [];
%     tmpxsp = superpixels(tmpx);
%     C = intersect(neighsp, tmpxsp);
%     if ~isempty(C)
%         elipos_x(i) = 0;        
%     end
% end
% elipos_x(elipos_x == 0) = [];

