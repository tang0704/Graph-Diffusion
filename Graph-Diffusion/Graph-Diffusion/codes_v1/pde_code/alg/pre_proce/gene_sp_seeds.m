function [seeds_sp] = gene_sp_seeds(seeds_location,superpixels)
% figure out which superpixel the 25 pixels belong to
% return index of superpixels, seeds_sp (25*1)

seeds_sp = zeros(length(seeds_location), 1);

for i = 1: length(seeds_location)
    x = int32(floor(seeds_location(i,1)));
    y = int32(floor(seeds_location(i,2)));
    seeds_sp(i) = superpixels(y, x); 
end    



