% function spAffinityMat = gene_weight( spAdjacentMat, sp_fea, theta,opts,neg_label_inds,pos_label_inds, n_ring)
function spAffinityMat = gene_weight( spAdjacentMat, sp_fea, sigmac, sigmas,opts,neg_label_inds,seeds_sp, n_ring, sp_center,Back_sp,Front_sp)
% get edges
% theta -- control the edge weight 
%  rgb_fea   
%  neg_label_inds   
%
% jjcao
%

if nargin < 7
    n_ring = 2;
end
spAffinityMat = spAdjacentMat;

%%
% %% compute affinity matrix structure: 
% for i = 2:n_ring % 2 or more rings of superpixels
%     spAffinityMat = spAffinityMat * spAdjacentMat;
% end
% if n_ring > 1
%     spAffinityMat(spAffinityMat>0) = 1;
%     spAffinityMat = spAffinityMat - eye(size(spAffinityMat));    
% end
% % attach addtional connection/edge
% spAffinityMat = append_connection_to_adjacent_matrix(spAffinityMat,neg_label_inds);
% spAffinityMat = append_connection_to_adjacent_matrix(spAffinityMat,pos_label_inds);

%% compute 1-ring & 2-ring of superpixels
sp_num = size(spAffinityMat,1);
edges=[];
for i=1:sp_num
    indext=[];
    ind=find(spAffinityMat(i,:)==1);
    
    % this for is for 2-ring
    if n_ring == 2
        for j=1:length(ind)
            indj=find(spAffinityMat(ind(j),:)==1);
            indext=[indext,indj];
        end
    end
    
    indext=[indext,ind];
    indext=indext((indext>i));
    indext=unique(indext);
    
    if(~isempty(indext))
        ed=ones(length(indext),2);
        ed(:,2)=i;
        ed(:,1)=indext;
        edges=[edges;ed];
    end
end



%% attach addtional connection/edge

edges1 = edges_between(neg_label_inds);
[meanDis,edges2] = edges_skeleton(seeds_sp, sp_center);
%  edges2=[];
edges = [edges; edges1; edges2];


Mid_sp=[];
if ~isempty(seeds_sp)
    for i=1:sp_num
        Mid_sp=[Mid_sp,i];
    end
    temp=[Back_sp', Front_sp'];
    Mid_sp=setdiff(Mid_sp, temp);
    Mid_sp=Mid_sp';
end


edges = remove_edges(edges, Back_sp,Front_sp);
% edges = remove_edges(edges, Back_sp,Mid_sp);


%% compute affinity matrix value
row = edges(:,1); col = edges(:,2);
ind = sub2ind(size(spAffinityMat), row, col); 

feat_dist_opt = getoptions(opts, 'feat_dist_opt', 'MVSSER');
if ~isstruct(sp_fea)
    feat_dist_opt = 'MVSSER';
end
switch feat_dist_opt
    case 'MVSSER'
        % feature
        tmpf = sum( (sp_fea(row,:) - sp_fea(col,:)).^2, 2);% feature
        tmpp = sum( (sp_center(row,:) - sp_center(col,:)).^2, 2);% distance
%         valDistances = sqrt(tmp)+eps;       
        tmpf = sqrt(tmpf)+eps;
        tmpp = sqrt(tmpp)+eps;
    case 'SDMVC'
        bin = fix((256-0.1)/sp_fea.binnum);
        sp_fea.R = (sp_fea.R+1).*bin;
        sp_fea.G = (sp_fea.G+1).*bin;
        sp_fea.B = (sp_fea.B+1).*bin;
        [L,a,b] = RGB2Lab(sp_fea.R,sp_fea.G,sp_fea.B);    
        tmp = patch_distance(sp_fea.pixelNumRGB,L,a,b,sp_fea.colorNumEachPatch, opts.npix, length(opts.npix));
        valDistances = tmp(ind)+eps;
    otherwise
        warning('not implement!');
end

% minVal = min(valDistances);
% valDistances=(valDistances-minVal)/(max(valDistances)-minVal);
mintmp1 = min(tmpf);
tmpf=(tmpf-mintmp1)/(max(tmpf)-mintmp1);    
mintmp2 = min(tmpp);
tmpp=(tmpp-mintmp2)/(max(tmpp)-mintmp2);
%     valDistances=normalize(valDistances); %Normalize to [0,1]. 



sigma=0.00255*meanDis;

if isempty(seeds_sp)
    sigma=0.125;
end    
% weights=exp(-theta*valDistances);%modified by Tang
% weights=exp(-tmpf/sigmas-tmpp/sigmac);
% weights=exp(-tmpf/sigmas);
weights=exp(-tmpf/sigmas-tmpp/sigma);


spAffinityMat(ind) = weights;
ind = sub2ind(size(spAffinityMat), col, row); 
spAffinityMat(ind) = weights;


if ~isempty(seeds_sp)
    for i=1:length(edges2)
        spAffinityMat(edges2(i,1),edges2(i,2))=1;
        spAffinityMat(edges2(i,2),edges2(i,1))=1;
    end    
end   


return;
end

function edges = edges_between(inds)
    %
    %    
    if isempty(inds)
        edges = [];
    end
    
    num = length(inds);   
    mat = tril(ones(num), -1);
    [row, col] = find(mat);
    edges = [inds(row), inds(col)];
    
end


function [meanDis, edges] = edges_skeleton(seeds_sp, sp_center)

%     (1, 2), (2, 21), (3, 21),
%     (4, 3), (5, 21), (6, 5), (7, 6), (8, 7), (9, 21),
%     (10, 9), (11, 10), (12, 11), (13, 1), (14, 13),
%     (15, 14), (16, 15), (17, 1), (18, 17), (19, 18),
%     (20, 19), (22, 23), (23, 8), (24, 25), (25, 12)

    edges=[];
    meanDis=1;
    if ~isempty(seeds_sp)
        skeleton=[1, 2; 2, 21; 3, 21;4, 3; 5, 21; 6, 5; 7, 6; 8, 7; 9, 21;
        10, 9; 11, 10; 12, 11; 13, 1; 14, 13;15, 14; 16, 15; 17, 1; 18, 17; 
         19, 18;20, 19; 22, 23; 23, 8; 24, 25; 25, 12];

        skeleton_num=size(skeleton,1);
        for i=1:skeleton_num
            x=skeleton(i,1);
            y=skeleton(i,2);
            a=seeds_sp(x);
            b=seeds_sp(y);
            edges=[edges;a, b];
        end   
        
        
        dis=0;
        for j=1:skeleton_num
            a=edges(j,1); b=edges(j,2);
            dis=dis+sqrt(sum((sp_center(a,:)-sp_center(b,:)).^2, 2));
        end    
        meanDis=dis/skeleton_num;
    end
  
end


function temp = remove_edges(edges, Back_sp,Front_sp)
    temp=edges;
    n=[];
    for i=1:length(edges)
       a=edges(i,1);b=edges(i,2);
       if ismember(a,Back_sp) 
           if ismember(b,Front_sp)
               %temp(i,:)=[];
               n=[n;i];
           end   
       elseif ismember(a,Front_sp) 
           if ismember(b,Back_sp)
%                temp(i,:)=[];
               n=[n;i];
           end 
       end    
        
    end  
   
    temp(n,:)=[];
end

function adjcMat = append_connection_to_adjacent_matrix(adjcMat,inds)
    %
    %    
    n = length(inds);    
    for i=1:n
        for j=i+1:n
            adjcMat(inds(i),inds(j))=1;
            adjcMat(inds(j),inds(i))=1;
        end
    end
end