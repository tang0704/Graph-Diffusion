function read_video(imdir, videoname)
% read videos, save videos as pic according to the frames
% return video names and the number of frames


video_num=length(videoname);

for j=1:1
%     current_viname=videoname(j).name;
    current_viname=videoname;
    video=VideoReader([imdir current_viname]); 
    frame_num=video.NumberOfFrames;
    frame_dir=[imdir 'img_' current_viname(1:end-4) '/'];
    
    
    if ~isdir(frame_dir)
        mkdir(frame_dir);
    end
    
    
    for i=1:frame_num
        frame=read(video,i);
        if mod(i,10)==1
            imwrite(frame, [frame_dir current_viname(1:end-4) '_' num2str(i) '.bmp'],'bmp');
        end    
        
    end 
end    

end