% Motion segmentation with Robust Shape Interaction Matrix Method with JBLD

clear; close all
addpath(genpath('../3rdParty'));
addpath(genpath('../matlab'));
dataName = 'Hopkins155';
dataPath = fullfile('~/research/data', dataName);
resPath = fullfile('../expData', 'res');
if ~exist(resPath, 'dir'), mkdir(resPath); end

func = {
    '[missrate,C1, label, gt] = SSCwrapper(X, s);';
    '[missrate, grp, CKSym, index] = ssc_JBLD(X, s);';
    '[missrate, grp, index] = SSC_View_JBLD_Obj(X, s, camID);';
    '[missrate, label, gt] = LRRwrapper(X, s);';
    '[missrate, label, gt] = LRR_JBLD(X, s);';
    '[missrate, grp, index] = LRR_View_JBLD_Obj(X, s, camID);';
    '[missrate, grp, bestRank, minNcutValue,W, index] = RSIM(X, s);';
    '[missrate, grp, bestRank, minNcutValue,W, index] = RSIM_JBLD(X, s, 4, 1);';
    '[missrate, grp, index] = RSIM_View_JBLD_Obj(X, s, 4, 1, camID);';
    '[missrate, grp, bestRank,W, index] = imprvRSIM(X, s, 4, 1, camID);';
    '[missrate, grp, bestRank,W, index] = imprvRSIM_JBLD(X, s, 4, 1, camID);';
    };

outputFile = 'output.txt';
fid = fopen(outputFile, 'w');

for funIndex = 6
% for funIndex = 1:length(func)    
    
    file = listFolder(dataPath);
    ii = 0;
    ii2 = 0;
    ii3 = 0;
    tic
    for i = 1:length(file)
%     for i = 82
        filePath = file{i};
        f = dir(fullfile(dataPath, filePath));
        foundValidData = false;
        for j = 1:length(f)
            if( ~isempty(strfind(f(j).name,'_truth.mat')) )
                ind = j;
                foundValidData = true;
                load(fullfile(dataPath, filePath, f(ind).name));
                if(max(s)==5)
                    foundValidData = false;
                end
                break;
            end
        end
        
        if(foundValidData)
            
            [s, ind] = sort(s); x = x(:, ind, :); y = y(:, ind, :);
            
%             % rotation and translation
%             theta = pi / 4; T = [300; 200]; ratio = 0.5;
%             [x,y,camID] = dataTransform(x, y, K, theta, T, ratio);
%             % rotation on half of the data
%             theta = pi / 4; ratio = 0.5;
%             [x, y, camID] = dataRotation(x, y, K, theta, ratio);
%             % translation on half of the data
%             T = [300; 200]; ratio = 0.5;
%             [x, y, camID] = dataTranslation(x, y, K, T, ratio);
            % delay on half of the data
            delay = 4; ratio = 0.5;
            [x, y, camID] = dataDelay(x, y, delay, ratio);
            
            [camID, idx] = sort(camID); x = x(:, idx, :); y = y(:, idx, :); s = s(idx);
            
%             % velocity
%             v = ones(size(x));
%             v(:,:,end) = [];
%             v(1:2,:,:) = diff(x(1:2,:,:), 1, 3);
%             x = v;
            
            N = size(x,2);
            F = size(x,3);
            D = 3*F;
            X = reshape(permute(x(1:3,:,:),[1 3 2]),D,N);	% note here the all-one rows are also included
            
            %             [missrate, grp, bestRank, minNcutValue,W, index] = RSIM(X, s);
            %             [missrate, grp, bestRank, minNcutValue,W, index] = RSIM_JBLD(X, s, 4, 1);
            %             [missrate, grp, bestRank, minNcutValue,W] = RSIM_JBLD2(X, s, 4, 1);
            %             [missrate, grp, bestRank, minNcutValue,W, index] = RSIM_JBLD3(X, s, 4, 1, camID);
            %             [missrate, grp, bestRank, minNcutValue,W, index] = RSIM_JBLD_oneshot(X, s, 4, 1, camID);
            %             [missrate, grp, bestRank, minNcutValue,W] = RSIM_JBLD4(X, s, 4, 1, camID);
            %             [missrate, grp, index] = RSIM_View_JBLD_Obj(X, s, 4, 1, camID);
            %             [missrate, grp, bestRank,W, index] = imprvRSIM(X, s, 4, 1, camID);
            %             [missrate, grp, bestRank, minNcutValue,W, index] = imprvRSIM_JBLD(X, s, 4, 1, camID);
            %             [missrate, grp, bestRank,W, index] = imprvRSIM_JBLD2(X, s, 4, 1, camID);
            %             [missrate, grp, bestRank,W, index] = imprvRSIM_JBLD_oneshot(X, s, 4, 1, camID);
            %             [missrate,C1, label, gt] = SSCwrapper(X, s);
            %             [missrate, grp, CKSym, index] = ssc_JBLD(X, s, camID);
            %             [missrate, grp, index] = SSC_View_JBLD_Obj(X, s, camID);
            
            eval(func{funIndex});
            
            
            if ~exist('label','var') && ~exist('gt','var')
%                 [~, labelIndex] = sort(index);
%                 label = grp * labelIndex';
%                 gt = s;
                nCluster = size(grp, 2);
                label = grp * (1:nCluster)';
                gt = index(s)';
            end
            m1 = nnz(label(camID==1)~=gt(camID==1));
            m2 = nnz(label(camID==2)~=gt(camID==2));
            tc = nnz(label~=gt);
            save(fullfile(resPath, ['label_', filePath]), 'label');
            clear label gt
            
            ii = ii+1;
            Missrate(ii) = missrate;
            M1(ii) = m1;
            M2(ii) = m2;
            TC(ii) = tc;
%                         disp([filePath ': ' num2str(100*Missrate(ii)) '%, dim:' num2str(bestRank) ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
            disp([filePath ': ' num2str(100*Missrate(ii)) '%, c1: ' ...
                num2str(M1(ii)) ', c2: ' num2str(M2(ii)) ', total error: ' num2str(TC(ii)) ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
%             disp([filePath ': ' num2str(100*Missrate(ii)) '%, nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
            if(max(s)==2)
                ii2 = ii2+1;
                Missrate2(ii2) = Missrate(ii);
            else
                ii3 = ii3+1;
                Missrate3(ii3) = Missrate(ii);
            end
        end
        %     end
    end
    time = toc;
    avgtime = time/ii;
    
    avgtol = mean(Missrate);
    medtol = median(Missrate);
    avgtwo = mean(Missrate2);
    medtwo = median(Missrate2);
    avgthree = mean(Missrate3);
    medthree = median(Missrate3);
    sumC1 = sum(M1);
    avgC1 = mean(M1);
    sumC2 = sum(M2);
    avgC2 = mean(M2);
    sumTC = sum(TC);
    avgTC = mean(TC);
    
    disp(['Results on ', dataName]);
    disp(['Function name: ', func{funIndex}]);
    disp(['Mean of two: ' num2str(100*avgtwo) '%' ', median of two: ' num2str(100*medtwo) '%;']);
    disp(['Mean of three: ' num2str(100*avgthree) '%' ', median of three: ' num2str(100*medthree) '%;']);
    disp(['Mean of all: ' num2str(100*avgtol) '%' ', median of all: ' num2str(100*medtol) '%;']);
    disp(['error # of cam1: ' num2str(sumC1) ', error # of cam2: ' num2str(sumC2) ', total error #: ' num2str(sumTC) '.']);
    disp(['avg error # of cam1: ' num2str(avgC1) ', avg error # of cam2: ' num2str(avgC2) ', avg total error #: ' num2str(avgTC) '.']);
    disp(['total running time: ', num2str(time), ' average running time per sequence: ', num2str(avgtime)]);
    disp([]);
    
    fprintf(fid, '%s\n', ['Results on ', dataName]);
    fprintf(fid, '%s\n', ['Function name: ', func{funIndex}]);
    fprintf(fid, '%s\n', ['Mean of two: ' num2str(100*avgtwo) '%' ', median of two: ' num2str(100*medtwo) '%;']);
    fprintf(fid, '%s\n', ['Mean of three: ' num2str(100*avgthree) '%' ', median of three: ' num2str(100*medthree) '%;']);
    fprintf(fid, '%s\n', ['Mean of all: ' num2str(100*avgtol) '%' ', median of all: ' num2str(100*medtol) '%;']);
    fprintf(fid, '%s\n', ['error # of cam1: ' num2str(sumC1) ', error # of cam2: ' num2str(sumC2) ', total error #: ' num2str(sumTC) '.']);
    fprintf(fid, '%s\n', ['avg error # of cam1: ' num2str(avgC1) ', avg error # of cam2: ' num2str(avgC2) ', avg total error #: ' num2str(avgTC) '.']);
    fprintf(fid, '%s\n', ['total running time: ', num2str(time), ' average running time per sequence: ', num2str(avgtime)]);
    fprintf(fid, '\n');
    
    res.index = funIndex;
    res.funName = func{funIndex};
    res.dataName = dataName;
    res.avgtwo = 100*avgtwo;
    res.medtwo = 100*medtwo;
    res.avgthree = 100*avgthree;
    res.medthree = 100*medthree;
    res.avgtol = 100*avgtol;
    res.medtol = 100*medtol;
    res.sumC1 = sumC1;
    res.sumC2 = sumC2;
    res.sumTC = sumTC;
    res.avgC1 = avgC1;
    res.avgC2 = avgC2;
    res.avgTC = avgTC;
    res.time = time;
    res.avgtime = avgtime;
    
    if ~exist('results', 'var')
        results = res;
    else
        results(end+1) = res;
    end
end
save results results
