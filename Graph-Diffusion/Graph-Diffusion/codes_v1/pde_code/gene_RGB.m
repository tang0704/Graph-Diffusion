function rgb = gene_RGB( s_num, spAdjcMat, seeds_sp, superpixels, input_im,sp_saliency)
% get RGB feature for every sp on one frame
% return final rgb(25,3)

adj_mat={s_num,1};
for i=1:s_num
    sp=seeds_sp(i);
    [ind]=find(spAdjcMat(sp,:)==1); 
    adj_mat{i,1}=ind;
end
 
rgb=zeros(s_num,3);
for i=1:s_num
    root_sp_ind=seeds_sp(i);
    group_sp=[];
    group_sp=[adj_mat{i,1} root_sp_ind];ֵ
    length_sp=length(group_sp);
    rgb_sp_r=0;rgb_sp_g=0;rgb_sp_b=0;
    
    for j=1:length_sp
        index_sp=group_sp(j);
        
        [row, col]=find(superpixels==index_sp);
        
        rgb_r=0;rgb_g=0;rgb_b=0;
        rgb_rr=zeros(length(row),1);
        rgb_gg=zeros(length(row),1);
        rgb_bb=zeros(length(row),1);
        for c=1:length(row)
            rgb_rr(c,1)=input_im(row(c),col(c),1);
            rgb_gg(c,1)=input_im(row(c),col(c),2);
            rgb_bb(c,1)=input_im(row(c),col(c),3);        
        end
        
        rgb_r=sum(rgb_rr);
        rgb_g=sum(rgb_gg);
        rgb_b=sum(rgb_bb);
        
        rgb_r=rgb_r/length(row);
        rgb_g=rgb_g/length(row);
        rgb_b=rgb_b/length(row);
        
        rgb_sp_r=rgb_sp_r+rgb_r*sp_saliency(index_sp);
        rgb_sp_g=rgb_sp_g+rgb_g*sp_saliency(index_sp);
        rgb_sp_b=rgb_sp_b+rgb_b*sp_saliency(index_sp);
    end
    
    
    rgb(i,1)=rgb_sp_r/length_sp;
    rgb(i,2)=rgb_sp_g/length_sp;
    rgb(i,3)=rgb_sp_b/length_sp;
end









