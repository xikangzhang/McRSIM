% Motion segmentation with Robust Shape Interaction Matrix Method with JBLD

clear; close all
addpath(genpath('../3rdParty'));
addpath(genpath('../matlab'));
dataPath = '~/research/data/Hopkins155';

rng('default');
file = listFolder(dataPath);
ii = 0;
ii2 = 0;
ii3 = 0;
% rng(3);
tic
for i = 1:length(file)
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
            
            for kk = 1:10
            [s, ind] = sort(s); x = x(:, ind, :);
            x3 = permute(x, [1 3 2]);
            N = size(x3, 3);
            F = size(x3, 2);
            D = 2 * F;
            X = reshape(x3(1:2,:,:), D, N);
            
            msk = (rand(F, N) > 0.05);
            indX = kron(msk, [1; 1]);
%             indX = kron(msk, [1; 1; 1]); indX(3:3:end, :) = 1;
            A = rand(D, N);
            X(indX==0) = A(indX==0);
            
%             [missrate,C1, label, gt] = SSCwrapper(X, s, 1);
            [missrate, grp, CKSym, index] = ssc_JBLD(X, s, 1);
%             [missrate, label, gt] = LRRwrapper(X, s, 1);
%             [missrate, label, gt] = LRR_JBLD(X, s);
            %             [missrate, grp, bestRank, minNcutValue,W, index] = RSIM(X, s);
%             [missrate, grp, bestRank, minNcutValue,W, index] = RSIM_JBLD(X, s, 4, 1);
            
            if ~exist('label','var') && ~exist('gt','var')
                nCluster = size(grp, 2);
                label = grp * (1:nCluster)';
                gt = index(s)';
            end
%             m1 = nnz(label(camID==1)~=gt(camID==1));
%             m2 = nnz(label(camID==2)~=gt(camID==2));
            tc = nnz(label~=gt);
            clear label gt
            
            ii = ii+1;
            Missrate(ii) = missrate;
%             M1(ii) = m1;
%             M2(ii) = m2;
            TC(ii) = tc;
%             disp([filePath ': ' num2str(100*Missrate(ii)) '%, dim:' num2str(bestRank) ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
%             disp([filePath ': ' num2str(100*Missrate(ii)) '%, c1: '  num2str(M1(ii)) ...
%                  'total error: ' num2str(TC(ii)) ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
%                 ', c2: ' num2str(M2(ii)) , ...
            disp([filePath ': ' num2str(100*Missrate(ii)) '%' ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
               
            if(max(s)==2)
                ii2 = ii2+1;
                Missrate2(ii2) = Missrate(ii);
            else
                ii3 = ii3+1;
                Missrate3(ii3) = Missrate(ii);
            end
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
% sumC1 = sum(M1);
% sumC2 = sum(M2);
sumTC = sum(TC);

disp('Results on Hopkins155 with gross corrupted entries')
disp(['Mean of two: ' num2str(100*avgtwo) '%' ', median of two: ' num2str(100*medtwo) '%;']);
disp(['Mean of three: ' num2str(100*avgthree) '%' ', median of three: ' num2str(100*medthree) '%;']);
disp(['Mean of all: ' num2str(100*avgtol) '%' ', median of all: ' num2str(100*medtol) '%;']);
% disp(['error # of cam1: ' num2str(sumC1) ', error # of cam2: ' num2str(sumC2) ', total error #: ' num2str(sumTC) '.']);
% disp(['error # of cam1: ' num2str(sumC1)  ', total error #: ' num2str(sumTC) '.']);
