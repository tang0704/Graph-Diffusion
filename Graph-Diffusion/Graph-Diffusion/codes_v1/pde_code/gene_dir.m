function [imdir, spdir, saldir,gdir, seeddir,rgbdir] = gene_dir(database,database_out)
% imdir -- test image path
% spdir -- the superpixel label file path
% saldir  -- the output path of the saliency map
% gdir -- the groundtruth path of image

% imdir=['../img/input/' database '/'];
% imdir=['../../test/'];
imdir=['E:/DATA/testRGB/tt/'];
% imdir=['C:/Users/tony/Desktop/error_video/aa/'];

if ~isdir(imdir)
    error('No image direction')
end

% spdir='../img/output/superpixels/';
spdir='../img/output/superpixel1/';
if ~isdir(spdir)
    mkdir(spdir);
end
 
% saldir=['../img/output/saliencymap/' database_out '/']; 
saldir=['../img/output/saliencyMap18/']; 
if ~isdir(saldir)
    mkdir(saldir);
end

gdir=['../img/G/' database '/']; 
if ~isdir(gdir)
    mkdir(gdir);
end

%skeleton����dir
% seeddir=['../../data1/']; 
% seeddir=['E:/DATA/skeleton/'];
seeddir=['E:/DATA/skeleton/'];
% seeddir=['C:/Users/Tang/Desktop/final/ttt/'];
% if ~isdir(seeddir)
%     mkdir(seeddir);
% end

rgbdir=['../../rgb/']; 
if ~isdir(rgbdir)
    mkdir(rgbdir);
end

