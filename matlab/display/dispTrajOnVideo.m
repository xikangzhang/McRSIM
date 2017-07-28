% display trajectories on videos

clear; close all; dbstop if error;
dataPath = '~/research/data/MultiViewMotionRaw';
seqPath = 'toy5305_seq1';
camName = 'cam2';
fileName = [seqPath, '_', camName '.avi'];
[~, f, ~] = fileparts(fileName);
load(fullfile(dataPath, seqPath, [f '_truth.mat']));

[s, ind] = sort(s); y = y(:, ind, :);
% load label;
% if strcmp(camName, 'cam1')
%     s = label(1:length(s));
% elseif strcmp(camName, 'cam2')
%     s = label(end-length(s)+1:end);
% end

% c = colormap('lines');
c = 'gym';

vidObj = VideoReader(fullfile(dataPath, seqPath, fileName));
currAxes = axes;
count = 1;
while hasFrame(vidObj)
    vidFrame = readFrame(vidObj);
    image(vidFrame, 'Parent', currAxes);
    currAxes.Visible = 'off';
    t = squeeze(y(1:2, :, count));
    hold on;
    for i = 1:max(s)
        plot(t(1, s==i), t(2, s==i), 'x', 'Color', c(i));
    end
    hold off;
    count = count + 1;
    pause(0.5);
end