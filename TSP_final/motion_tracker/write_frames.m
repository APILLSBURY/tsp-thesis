function frame_rate = write_frames(movie_name, output_folder, img_suffix)
    disp('Writing frames to disk...');
    videoObject = VideoReader(movie_name);
    frame_rate = videoObject.FrameRate;
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    if numel(dir([output_folder '*.' img_suffix]))==0
        f=1;
        while hasFrame(videoObject)
            fprintf('Writing frame %d to disk\n', f);
            thisFrame = readFrame(videoObject);
            frame_name = fullfile(output_folder, sprintf(['frame%4.4d.' img_suffix], f));
            imwrite(thisFrame, frame_name, img_suffix);
            f=f+1;
        end
    end
end
