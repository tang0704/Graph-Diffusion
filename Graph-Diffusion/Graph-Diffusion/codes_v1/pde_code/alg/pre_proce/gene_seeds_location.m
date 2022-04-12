function [seeds_location1, seeds_location2, real_p_num, numerror,numless1,data0] = gene_seeds_location(fid,index,whole_frame_num,v_name)
% read location data of 25 seeds
% return RGB location, seeds_location (25*2)
% per_num -- number of people in this frame
% seeds_location1 -- location of the first person
% seeds_location2 -- location of the second person(if exsit)
% index -- the index(th) frame
% whole_frame_num -- the number of the whole frames
% numerror -- if the number of people in pic >2 , numeeor=1
% real_p_num -- real number of people in one frame
% numless1 -- real number is 2 while the number in data is 1.it's a flag

next = index*10+1;


seeds_location1 = zeros(25,2);
seeds_location2 = zeros(25,2);

str=fgetl(fid);
per_num=str2num(str);


  
real_p_num=0;
numless1=0;
numerror=0;

data0=0;
if per_num==0
   data0=1; 
   return;
end  

if per_num >2
    numerror=1;
    return;
end    


real_flag_num=str2num(v_name(18:20));
if real_flag_num<=49
    real_p_num=1;
elseif real_flag_num>49 &&real_flag_num<=60
    real_p_num=2;
end    

if real_p_num==2 && per_num==1
   numless1=1;
   return;
end    



fgetl(fid);fgetl(fid);
temp1 = zeros(1, 12);
temp2 = zeros(1, 12);

for j = 1:25
    line = fgetl(fid);
    part = regexp(line, '\s+');
    content = zeros(1, 12);
    
    for i = 1:length(part)+1
        
        if i == 1
            content(1, i)= str2num(line(1:part(i)));
        elseif i == length(part)+1
            content(1, i)= str2num(line(part(i-1)+1:end));
        else
            start=part(i-1)+1; endd=part(i)-1;
            content(1, i)= str2num(line(start:endd));
        end
    end
    
    if content(1,6) > 1920
        seeds_location1(j, 1) = 1920;
    elseif isnan(content(1,6))
        seeds_location1(j, 1) = temp1(1,6);
    else
        seeds_location1(j, 1) = content(1, 6);
    end    
    
    if content(1,7) > 1080
        seeds_location1(j, 2)=1080;
    elseif isnan(content(1,6))
        seeds_location1(j, 1) = temp1(1,7);
    else
        seeds_location1(j, 2) = content(1, 7);
    end    
    temp1 = content;
    
end



if per_num==2 && real_p_num==1
    fgetl(fid);fgetl(fid);
    ii=25;
    while(ii>0)
        fgetl(fid);
        ii=ii-1;
    end
end




if real_p_num>1
    fgetl(fid);fgetl(fid);
    for j = 1:25
        line = fgetl(fid);
        part = regexp(line, '\s+');
        content = zeros(1, 12);
        
        for i = 1:length(part)+1
            
            if i == 1
                content(1, i)= str2num(line(1:part(i)));
            elseif i == length(part)+1
                content(1, i)= str2num(line(part(i-1)+1:end));
            else
                start=part(i-1)+1; endd=part(i)-1;
                content(1, i)= str2num(line(start:endd));
            end
        end

        if content(1,6) > 1920
            seeds_location2(j, 1) = 1920;
        elseif isnan(content(1,6))
            seeds_location2(j, 1) = temp2(1,6);
        else
            seeds_location2(j, 1) = content(1, 6);
        end
        
        if content(1,7) > 1080
            seeds_location2(j, 2)=1080;
        elseif isnan(content(1,6))
            seeds_location2(j, 1) = temp2(1,7);
        else
            seeds_location2(j, 2) = content(1, 7);
        end
        temp2 = content;
    end
end


if next<whole_frame_num
    for i=1:9
        p_num=fgetl(fid);
        
        if p_num==0
           data=0;
           return;
        end    
        
        if real_p_num==2 && p_num==1
            numless1=1;
            return;
        end
        
        if isa(p_num, 'double')
           p_num=num2str(per_num);
        end    
%         fgetl(fid);
%         p_num=1;

%         if str2num(p_num)==1
        if real_p_num==1
            a=fgetl(fid);fgetl(fid);
            ii=25;
            while(ii>0)
                fgetl(fid);
                ii=ii-1;
            end
            
            
            if str2num(p_num)==2
                fgetl(fid);fgetl(fid);
                ii=25;
                while(ii>0)
                    fgetl(fid);
                    ii=ii-1;
                end
            end
            
%             fgetl(fid);fgetl(fid);
%             ii=25;
%             while(ii>0)
%                 fgetl(fid);
%                 ii=ii-1;
%             end
            
            
            
            
        else
            fgetl(fid);fgetl(fid);
            ii=52;
            while(ii>0)
                fgetl(fid);
                ii=ii-1;
            end
        end
    end
end    


