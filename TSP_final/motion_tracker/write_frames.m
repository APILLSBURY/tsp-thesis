function frame_rate = write_frames(movie_name, output_folder, img_suffix)
    disp('Writing frames to disk...');
    videoObject = VideoReader(movie_name);
    frame_rate = videoObject.FrameRate;
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    if numel(dir([output_folder '*.' img_suffix]))==0
        frames = read(videoObject);
        for f = 1:size(frames, 4)
            fprintf('Writing frame %d to disk\n', f);
            frame_name = fullfile(output_folder, sprintf(['frame%4.4d.' img_suffix], f));
            imwrite(frames(:,:,:,f), frame_name, img_suffix);
        end
    end
end
