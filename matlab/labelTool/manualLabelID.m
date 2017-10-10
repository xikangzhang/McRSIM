function manualLabelID

close all;
dataPath = '~/research/data/RSL60';
seqPath = 'toy2_038';
camName = 'cam2';
fileName = [seqPath, '_', camName, '.avi'];
[~, f, ~] = fileparts(fullfile(dataPath, seqPath, fileName));
load(fullfile('../expData', [f '_xyK.mat']));
loc = y(1:2, :, 1)';

vidObj = VideoReader(fullfile(dataPath, seqPath, fileName));
vidFrame = readFrame(vidObj);
displayPoints(vidFrame, loc);

% outlier removal
fprintf('Now use mouse to remove the outliers.\n');
r = [];
while isempty(r)
    rect = getrect;
    bb = [rect(1), rect(1)+rect(3)-1, rect(2), rect(2)+rect(4)-1];
    isInBbox = (loc(:, 1) >= bb(1) & loc(:, 1) <= bb(2) & loc(:, 2) >= bb(3) & loc(:, 2) <= bb(4));
    loc(isInBbox, :) = [];
    y(:, isInBbox, :) = [];
    x(:, isInBbox, :) = [];
    displayPoints(vidFrame, loc);
    r = input('Press Enter to continue. Press q if you are done. \n','s');
end

% initialize ground truth label s
nPoint = size(y, 2);
s = zeros(nPoint, 1);
index = (1:nPoint)';

m = input('Input number of motions:  ');
if ~isnumeric(m)
    error('Input must be a number.\n');
end

for i = 1:m
    fprintf('Now use mouse to select the points in Cluster %d. \n', i);
    r = [];
    while isempty(r)
        rect = getrect;
        bb = [rect(1), rect(1)+rect(3)-1, rect(2), rect(2)+rect(4)-1];
        isInBbox = (loc(:, 1) >= bb(1) & loc(:, 1) <= bb(2) & loc(:, 2) >= bb(3) & loc(:, 2) <= bb(4));
        s(index(isInBbox)) = i;
        loc(isInBbox, :) = [];
        index(isInBbox) = [];
        displayPoints(vidFrame, loc);
        r = input('Press Enter to continue. Press q if you are done. \n','s');
    end
    
end

save(fullfile('../expData', [f '_truth']), 'x', 'y', 'K', 's');

end

function displayPoints(vidFrame, loc)

imshow(vidFrame);
hold on;
title('Detected features');
plot(loc(:, 1), loc(:, 2), 'g+');
hold off;

end