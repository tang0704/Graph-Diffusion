function [cntSeg, npixSeg, sp_inds, contourSeg]= compute_sp_info(spImg, nSeg)
%
% center coordinates & number of pixel of each segment
%
% jjcao @ 2012

cntSeg = zeros(nSeg,2);
npixSeg = zeros(nSeg,1);
sp_inds=cell(nSeg,1);
[m, n] = size(spImg);

if nargout > 3
    imsize = size(spImg);
    contourSeg = cell(nSeg,1);
end
        
%%
for j=1:nSeg,
    sp_inds{j}=find(spImg==j);
    [r c] = ind2sub([m,n], sp_inds{j});
        
    npixSeg(j) = length(r);
    cntSeg(j,:) = round([mean(r) mean(c)]);
    
    if nargout > 3 
        ind = sub2ind(imsize, r, c);
        BW = zeros(imsize);
        BW(ind) = 1;
        B = bwboundaries(BW,'noholes');
        contourSeg{j} = B{1};

    %     imshow(label2rgb(spImg, @jet, [.5 .5 .5]));  hold on;
    %     for k = 1:length(B)
    %         boundary = B{k};
    %         plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
    %     end
    end
end
    
end