% try KLT tracker in CV toolbox

% display trajectories on videos

clear; close all; dbstop if error;
dataPath = '~/research/data/MultiViewMotionRaw';
seqPath = 'toy5305_seq1';
camName = 'cam1';
fileName = [seqPath, '_', camName, '.avi'];

% [~,f,~] = fileparts(fileName);
% load(fullfile('../expData', [f '_xyK.mat']));
% points = squeeze(y(1:2, :, 1))';



vidObj = VideoReader(fullfile(dataPath, seqPath, fileName));
vidFrame = readFrame(vidObj);
wid = size(vidFrame, 2);
hgt = size(vidFrame, 1);

points = detectMinEigenFeatures(rgb2gray(vidFrame), 'MinQuality', 0.01);
% points = detectSURFFeatures(rgb2gray(vidFrame));
points = points.Location;
points = manualRefine(vidFrame, points);
% imshow(vidFrame), hold on, title('Detected features');
% plot(points);

pointTracker = vision.PointTracker('MaxBidirectionalError', 2);
initialize(pointTracker, points, vidFrame);

d = 2;
nFrame = 100;
nPoint = size(points, 1);
tracks = zeros(d, nPoint, nFrame);
tracks(:, :, 1) = points.';

figure;
currAxes = axes;
count = 2;
oldPoints = points;
while hasFrame(vidObj)
    
    vidFrame = readFrame(vidObj);
    image(vidFrame, 'Parent', currAxes);
    currAxes.Visible = 'off';

    [points, isFound] = step(pointTracker, vidFrame);
    visiblePoints = points(isFound, :);
    oldInliers = oldPoints(isFound, :);
    
    tracks(:, ~isFound, :) = []; % remove lost tracks
    tracks(:, :, count) =  visiblePoints.';
    
    if size(visiblePoints, 1) >= 2 % need at least 2 points
        oldPoints = visiblePoints;
        setPoints(pointTracker, oldPoints);
    end
    hold on;plot(visiblePoints(:, 1), visiblePoints(:, 2), 'gx');hold off;
    
    count = count + 1;
    count
    pause(0.1);
end
tracks(:, :, count:end) = [];

% Clean up
release(pointTracker);

y = cat(1, tracks, ones(1, size(tracks, 2), size(tracks, 3)));
K = [wid/2, 0 wid/2; 0, hgt/2, hgt/2; 0, 0, 1];
x = y;
for i = 1:size(y, 3)
    x(:, :, i) = K \ squeeze(y(:, :, i));
end
[~,f,~] = fileparts(fileName);
save(fullfile('../expData', [f '_xyK.mat']), 'x', 'y', 'K');