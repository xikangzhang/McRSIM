% load label, save videos with labels superimposed

% display trajectories on videos

clear; close all; dbstop if error;
dataPath = '~/research/data/MultiViewMotion2';
% dataPath = '~/research/data/Hopkins155';
labelPath = '../expData/res';
camList = {'cam1', 'cam2'};
% camList = {''};
videoOutPath = '../expData/outputVideos';
if ~exist(videoOutPath, 'dir'), mkdir(videoOutPath); end

seqPathList = listFolder(dataPath);
for i = 1:length(seqPathList)
% for i = 82
    seqPath = seqPathList{i};
    load(fullfile(labelPath, ['label_' seqPath]), 'label');
    for j = 1:length(camList)
        camName = camList{j};
        fileName = [seqPath, '_', camName '.avi'];
%         fileName = [seqPath, '.avi'];
        [~, f, ~] = fileparts(fileName);
        load(fullfile(dataPath, seqPath, [f '_truth.mat']));
        [s, ind] = sort(s); y = y(:, ind, :);
        
        if strcmp(camName, 'cam1') || isempty(camName)
            s = label(1:length(s));
        elseif strcmp(camName, 'cam2')
            s = label(end-length(s)+1:end);
        end
        
        % c = colormap('lines');
        c = 'gym';
%         c = 'bbb';
        
        vidIn = VideoReader(fullfile(dataPath, seqPath, fileName));
        vidOut = VideoWriter(fullfile(videoOutPath, [f, '_out.avi']));
        open(vidOut);
        currAxes = axes;
        count = 1;
        while hasFrame(vidIn)
            vidFrame = readFrame(vidIn);
            image(vidFrame, 'Parent', currAxes);
            currAxes.Visible = 'off';
            if count > size(y, 3), continue; end
            t = squeeze(y(1:2, :, count));
            hold on;
            for k = 1:max(s)
                plot(t(1, s==k), t(2, s==k), 'x', 'Color', c(k));
%                 plot(t(1, s~=label), t(2, s~=label), 'x', 'Color', c(k));
            end
            hold off;
            currFrame = getframe;
            writeVideo(vidOut, currFrame);
            count = count + 1;
            pause(0.01);
        end
        close(vidOut);
    end
end