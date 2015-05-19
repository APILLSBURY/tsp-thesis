function write_video(ims, vid_name, frame_rate)
    disp('Writing video...');
    % Create a VideoWriter object to write the video out to a new, different file.
	writerObj = VideoWriter(vid_name);
    writerObj.FrameRate = frame_rate;
	open(writerObj);
	for frame = 1 : size(ims, 4)
        fprintf('Writing frame %d/%d\n', frame, size(ims, 4));
		writeVideo(writerObj, ims(:,:,:,frame));
	end
	close(writerObj);
    disp('Done writing video!');
end