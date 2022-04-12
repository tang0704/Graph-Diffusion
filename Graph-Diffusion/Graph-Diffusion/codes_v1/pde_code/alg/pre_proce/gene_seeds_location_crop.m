function [input_im,seeds_location] = gene_seeds_location_crop(input_im,seeds_location)

% return 25 seeds location after crop

[borderx,bordery,c]=size(input_im);


minx=min(seeds_location(:,1));
miny=min(seeds_location(:,2));

maxx=max(seeds_location(:,1));
maxy=max(seeds_location(:,2));

lengthx=maxx-minx;
lengthy=maxy-miny;

alpha=0.35;
belta=0.24;

minx=minx-lengthx*alpha;
miny=miny-lengthy*belta;

if minx<0
   minx=0; 
end    
if miny<0
   miny=0; 
end  

lengthx=lengthx+lengthx*alpha*2;
lengthy=lengthy+lengthy*belta*2;
% maxx=maxx+lengthx*alpha*2;
% maxy=maxy+lengthy*alpha*2;


for i=1:length(seeds_location)
   seeds_location(i, 1)=seeds_location(i, 1)- minx+1;%x
   %·ÀÖ¹Ô½½ç
   if seeds_location(i, 1)>bordery
      seeds_location(i, 1)=bordery; 
   elseif seeds_location(i, 1)<1
       seeds_location(i, 1)=1;
   end    
   seeds_location(i, 2)=seeds_location(i, 2)- miny+1;%y
   if seeds_location(i, 2)>borderx
      seeds_location(i, 2)=borderx; 
   elseif seeds_location(i, 2)<1
       seeds_location(i, 2)=1;
   end 
end    

% figure('name', '11111222222222');
input_im=imcrop(input_im, [minx,miny,lengthx,lengthy]);
% imshow(input_im);


