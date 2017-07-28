% Motion segmentation with Robust Shape Interaction Matrix Method with JBLD

clear; close all
addpath(genpath('../3rdParty'));
addpath(genpath('../matlab'));
dataPath = '~/research/data/MultiViewMotion2';


file = listFolder(dataPath);
ii = 0;
ii2 = 0;
ii3 = 0;
tic
for i = 1:length(file)
% for i = 14
    [~, f, ~] = fileparts(file{i});
    load(fullfile(dataPath, file{i}, [f '_cam1_truth.mat']));
    x1 = x; s1 = s;
    [s1, ind] = sort(s1); x1 = x1(:, ind, :);
    load(fullfile(dataPath, file{i}, [f '_cam2_truth.mat']));
    x2 = x; s2 = s;
    [s2, ind] = sort(s2); x2 = x2(:, ind, :);
    x = cat(2, x1, x2); s = [s1; s2]; camID = [ones(length(s1), 1); 2*ones(length(s2), 1)];
    
%     x = x2; s = s2;
%     theta =  pi/4;
%     temp = x;
%     for k = 1:size(x, 3)
%         temp(1:2, :, k) = rotation(x(1:2, :, k), theta);
%     end
% %     temp = x([2 1 3], :, :);
%     z = cat(2, x, temp); x = z; camID = [ones(length(s), 1); 2*ones(length(s), 1)]; s = [s; s];
    
%     [s, ind] = sort(s); x = x(:, ind, :);

            N = size(x,2);
            F = size(x,3);
            D = 3*F;
            
            X = reshape(permute(x(1:3,:,:),[1 3 2]),D,N);	% note here the all-one rows are also included
            
%             [missrate, grp, bestRank, minNcutValue,W, index] = RSIM(X, s);
%             [missrate, grp, bestRank, minNcutValue,W, index] = RSIM_JBLD(X, s, 4, 1);
%             [missrate, grp, bestRank, minNcutValue,W] = RSIM_JBLD2(X, s, 4, 1);
%             [missrate, grp, bestRank, minNcutValue,W, index] = RSIM_JBLD3(X, s, 4, 1, camID);
%             [missrate, grp, bestRank, minNcutValue,W, index] = RSIM_JBLD_oneshot(X, s, 4, 1, camID);
%             [missrate, grp, bestRank, minNcutValue,W1, index] = RSIM_JBLD_oneshot2(X, s, 4, 1, camID);
%             [missrate, grp, bestRank, minNcutValue,W] = RSIM_JBLD4(X, s, 4, 1, camID);
            [missrate, grp, index] = RSIM_View_JBLD_Obj(X, s, 4, 1, camID);
%             [missrate, grp, bestRank,W, index] = imprvRSIM_JBLD(X, s, 4, 1, camID);
%             [missrate, grp, bestRank,W, index] = imprvRSIM_JBLD2(X, s, 4, 1, camID);
%             [missrate, grp, bestRank,W, index] = imprvRSIM_JBLD_oneshot(X, s, 4, 1, camID);

            nCluster = size(grp, 2);
            label = grp * (1:nCluster)';
            gt = index(s)';
            m1 = nnz(label(camID==1)~=gt(camID==1));
            m2 = nnz(label(camID==2)~=gt(camID==2));
            tc = nnz(label~=gt);
            
            ii = ii+1;
            Missrate(ii) = missrate;
            M1(ii) = m1;
            M2(ii) = m2;
            TC(ii) = tc;
%             disp([filePath ': ' num2str(100*Missrate(ii)) '%, dim:' num2str(bestRank) ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
            disp([file{i} ': ' num2str(100*Missrate(ii)) '%, c1: ' ...
                num2str(M1(ii)) ', c2: ' num2str(M2(ii)) ', total error: ' num2str(TC(ii)) ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
            if(max(s)==2)
                ii2 = ii2+1;
                Missrate2(ii2) = Missrate(ii);
            else
                ii3 = ii3+1;
                Missrate3(ii3) = Missrate(ii);
            end
end
time = toc;
avgtime = time/ii

avgtol = mean(Missrate);
medtol = median(Missrate);
avgtwo = mean(Missrate2);
medtwo = median(Missrate2);
avgthree = mean(Missrate3);
medthree = median(Missrate3);
sumC1 = sum(M1);
sumC2 = sum(M2);
sumTC = sum(TC);

disp('Results on MultiViewMotion')
disp(['Mean of all: ' num2str(100*avgtol) '%' ', median of all: ' num2str(100*medtol) '%;']);
disp(['Mean of two: ' num2str(100*avgtwo) '%' ', median of two: ' num2str(100*medtwo) '%;']);
disp(['Mean of three: ' num2str(100*avgthree) '%' ', median of three: ' num2str(100*medthree) '%;']);
disp(['error # of cam1: ' num2str(sumC1) ', error # of cam2: ' num2str(sumC2) ', total error #: ' num2str(sumTC) '.']);