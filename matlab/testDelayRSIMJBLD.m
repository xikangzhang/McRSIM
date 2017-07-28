% Motion segmentation with Robust Shape Interaction Matrix Method with JBLD

clear; close all
addpath(genpath('../3rdParty'));
addpath(genpath('../matlab'));
dataPath = '~/research/data/Hopkins155';

% file = dir(dataPath);
file = listFolder(dataPath);
ii = 0;
ii2 = 0;
ii3 = 0;
tic
for i = 1:length(file)
% for i = 12
%     if( (file(i).isdir == 1) && ~strcmp(file(i).name,'.') && ~strcmp(file(i).name,'..') )
%         filePath = file(i).name;
        filePath = file{i};
        % 		eval(['cd ' filePath]);
        f = dir(fullfile(dataPath, filePath));
        foundValidData = false;
        for j = 1:length(f)
            if( ~isempty(strfind(f(j).name,'_truth.mat')) )
                ind = j;
                foundValidData = true;
                % 				eval(['load ' f(ind).name]);
                load(fullfile(dataPath, filePath, f(ind).name));
                if(max(s)==5)
                    foundValidData = false;
                end
                break;
            end
        end
%         cd ..
        
        if(foundValidData)
            
            
            [s, ind] = sort(s);
            x = x(:, ind, :);
%             % x2 is delay of x1
%             nFrame = size(x, 3);
%             delay = 4;
%             x1 = x(:,:,1:nFrame-delay);
%             x2 = x(:,:,delay+1:nFrame);
% %             z = cat(2, x1, x2); x = z; camID = [ones(length(s), 1); 2*ones(length(s), 1)]; s = [s; s];
%             x = x1;
%             
%             N = size(x,2);
%             F = size(x,3);
%             D = 3*F;
%             
%             X = reshape(permute(x(1:3,:,:),[1 3 2]),D,N);	% note here the all-one rows are also included
            
            delay =  2;
            x3 = permute(x, [1 3 2]);
            N = size(x3, 3);
            F = size(x3, 2);
            D = 3 * (F-delay);
            
            rng(0);
            ind = (rand(1, N) > 0.5);
            camID = zeros(N, 1);
            x1 = x3(:, 1:F-delay, ind); camID(ind) = 1;
            x2 = x3(:, delay+1:F, ~ind); camID(~ind) = 2;
            x3 = cat(3, x1, x2); s = [s(ind); s(~ind)]; camID = [camID(ind); camID(~ind)];
            X = reshape(x3,D,N);
            
%             D = size(x, 1);
%             N = size(x, 2);
%             F = size(x, 3);
%             x4 = zeros(D, N, F-delay);
%             rng(0);
%             ind = (rand(1, N) > 0.5);
%             camID = zeros(N, 1);
%             x4(:, ind, :) = x(:, ind, 1:F-delay); camID(ind) = 1;
%             x4(:, ~ind, :) = x(:, ~ind, delay+1:F); camID(~ind) = 2;
%             x = x4;
%             [camID, idx] = sort(camID); x = x(:, idx, :); s = s(idx);
% %             x = cat(2, x1, x2); s = [s(ind); s(~ind)]; camID = [camID(ind); camID(~ind)];
% %             X = reshape(x3,D,N);
%             
%             N = size(x,2);
%             F = size(x,3);
%             D = 3*F;
%             
%             X2 = reshape(permute(x(1:3,:,:),[1 3 2]),D,N);	% note here the all-one rows are also included

%             [missrate, grp, bestRank, minNcutValue,W, index] = RSIM(X, s);
%             [missrate, grp, bestRank, minNcutValue,W, index] = RSIM_JBLD(X, s, 4, 1);
%             [missrate, grp, bestRank, minNcutValue,W] = RSIM_JBLD2(X, s, 4, 1);
%             [missrate, grp, bestRank, minNcutValue,W, index] = RSIM_JBLD3(X, s, 4, 1, camID);
%             [missrate, grp, bestRank, minNcutValue,W] = RSIM_JBLD4(X, s, 4, 1, camID);
%             [missrate, grp, index] = RSIM_View_JBLD_Obj(X, s, 4, 1, camID);
            [missrate, grp, bestRank, minNcutValue,W, index] = RSIM_JBLD_oneshot(X, s, 4, 1, camID);
%             [missrate, grp, bestRank,W, index] = imprvRSIM_JBLD2(X, s, 4, 1, camID);
            
%             x3 = x3(1:2,:,:); D = 2 * size(x3, 2); X = reshape(x3, D, N);
%             r = 0; affine = true; alpha = 800; outlier = false; rho = 0.7;
%             [missrate,C1, label, gt] = SSC(X,r,affine,alpha,outlier,rho,s);
            
            if ~exist('label','var') && ~exist('gt','var')
                nCluster = size(grp, 2);
                label = grp * (1:nCluster)';
                gt = index(s)';
            end
            m1 = nnz(label(camID==1)~=gt(camID==1));
            m2 = nnz(label(camID==2)~=gt(camID==2));
            tc = nnz(label~=gt);
            clear label gt;
            
            ii = ii+1;
            Missrate(ii) = missrate;
            M1(ii) = m1;
            M2(ii) = m2;
            TC(ii) = tc;
%             disp([filePath ': ' num2str(100*Missrate(ii)) '%, dim:' num2str(bestRank) ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
            disp([filePath ': ' num2str(100*Missrate(ii)) '%, c1: ' ...
                num2str(M1(ii)) ', c2: ' num2str(M2(ii)) ', total error: ' num2str(TC(ii)) ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
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

disp('Results on Hopkins155')
disp(['Mean of all: ' num2str(100*avgtol) '%' ', median of all: ' num2str(100*medtol) '%;']);
disp(['Mean of two: ' num2str(100*avgtwo) '%' ', median of two: ' num2str(100*medtwo) '%;']);
disp(['Mean of three: ' num2str(100*avgthree) '%' ', median of three: ' num2str(100*medthree) '%.']);
disp(['error # of cam1: ' num2str(sumC1) ', error # of cam2: ' num2str(sumC2) ', total error #: ' num2str(sumTC) '.']);