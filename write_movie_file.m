% combine still images into movie

% If your image file starts with image_1.png, image_2.png and so on ...
% and live in the current folder.
folder = pwd; % Or wherever you want.
video = VideoWriter('testvideo.mp4'); % Create the video object.
video.FrameRate = 1;
open(video); % Open the file for writing
N=8; % Where N is the separate number of PNG image files.
for k = 1 : N 
    I = imread(fullfile(folder, 'volume', 'plots', sprintf('area_ic5_sm1000_designmat1_SAgeSexScanSoft_age_%g.0.png', k+8))); % Read the next image from disk.
    writeVideo(video,I); % Write the image to file.
end
close(video); 

% /home/d9smith/projects/age/volume/plots/area_ic5_sm1000_designmat1_SAgeSexScanSoft_age_9.0.png