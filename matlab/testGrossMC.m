% Motion segmentation on Hopkins 155 (Half trajectories rotated) with gross contamination

clear; close all
addpath(genpath('../3rdParty'));
addpath(genpath('../matlab'));
dataPath = '~/research/data/Hopkins155';

func = {
    '[missrate, grp, bestRank, minNcutValue,W, index] = RSIM(X, s);';
    '[missrate, grp, bestRank, minNcutValue,W, index] = RSIM_MDD(X, s, 4, 1, 1);';
    '[missrate, grp, bestRank,W, index] = McRSIM(X, s, 4, 1, camID);';
    '[missrate, grp, bestRank,W, index] = McRSIM_MDD(X, s, 4, 1, camID, 1);';
    };

outputFile = 'outputGrossMC.txt';
fid = fopen(outputFile, 'w');

for funIndex = 1:length(func)
    rng('default');
    file = listFolder(dataPath);
    ii = 0;
    ii2 = 0;
    ii3 = 0;
    tic
    for i = 1:length(file)
        for kk = 1:10
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
                
                % randomly pick half the points and rotate theta
                [s, ind] = sort(s); x = x(:, ind, :); y = y(:, ind, :);
                theta = pi / 4; ratio = 0.5;
                [x, y, camID] = dataRotation(x, y, K, theta, ratio);
                [camID, idx] = sort(camID); x = x(:, idx, :); y = y(:, idx, :); s = s(idx);
                
                N = size(x,2);
                F = size(x,3);
                D = 2*F;
                
                X = reshape(permute(x(1:2,:,:),[1 3 2]),D,N);
                
                rng(kk);
                msk = (rand(F, N) > 0.05);
                indX = kron(msk, [1; 1]);
                A = rand(D, N);
                X(indX==0) = A(indX==0);
                
                eval(func{funIndex});
                
                if ~exist('label','var') && ~exist('gt','var')
                    nCluster = size(grp, 2);
                    label = grp * (1:nCluster)';
                    gt = index(s)';
                end
                m1 = nnz(label(camID==1)~=gt(camID==1));
                m2 = nnz(label(camID==2)~=gt(camID==2));
                tc = nnz(label~=gt);
                clear label gt x y
                
                ii = ii+1;
                Missrate(ii) = missrate;
                M1(ii) = m1;
                M2(ii) = m2;
                TC(ii) = tc;
                disp([filePath ': ' num2str(100*Missrate(ii)) '%, dim:' num2str(bestRank) ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
                %             disp([filePath ': ' num2str(100*Missrate(ii)) '%, c1: '  num2str(M1(ii)) ...
                %                  'total error: ' num2str(TC(ii)) ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
                %                 ', c2: ' num2str(M2(ii)) , ...
                
                if(max(s)==2)
                    ii2 = ii2+1;
                    Missrate2(ii2) = Missrate(ii);
                else
                    ii3 = ii3+1;
                    Missrate3(ii3) = Missrate(ii);
                end
                
            end
        end
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
    
    disp('Results on Hopkins155, half rotated 45 degree, with gross corrupted entries');
    disp(['Function name: ', func{funIndex}]);
    disp(['Mean of two: ' num2str(100*avgtwo) '%' ', median of two: ' num2str(100*medtwo) '%;']);
    disp(['Mean of three: ' num2str(100*avgthree) '%' ', median of three: ' num2str(100*medthree) '%;']);
    disp(['Mean of all: ' num2str(100*avgtol) '%' ', median of all: ' num2str(100*medtol) '%;']);
    disp(['error # of cam1: ' num2str(sumC1) ', error # of cam2: ' num2str(sumC2) ', total error #: ' num2str(sumTC) '.']);
    disp(['error # of cam1: ' num2str(sumC1)  ', total error #: ' num2str(sumTC) '.']);
    
    fprintf(fid, '%s\n', 'Results on Hopkins155, half rotated 45 degree, with gross corrupted entries');
    fprintf(fid, '%s\n', ['Function name: ', func{funIndex}]);
    fprintf(fid, '%s\n', ['Mean of two: ' num2str(100*avgtwo) '%' ', median of two: ' num2str(100*medtwo) '%;']);
    fprintf(fid, '%s\n', ['Mean of three: ' num2str(100*avgthree) '%' ', median of three: ' num2str(100*medthree) '%;']);
    fprintf(fid, '%s\n', ['Mean of all: ' num2str(100*avgtol) '%' ', median of all: ' num2str(100*medtol) '%;']);
    fprintf(fid, '%s\n', ['error # of cam1: ' num2str(sumC1) ', error # of cam2: ' num2str(sumC2) ', total error #: ' num2str(sumTC) '.']);
    fprintf(fid, '%s\n', ['avg error # of cam1: ' num2str(avgC1) ', avg error # of cam2: ' num2str(avgC2) ', avg total error #: ' num2str(avgTC) '.']);
    fprintf(fid, '%s\n', ['total running time: ', num2str(time), ', average running time per sequence: ', num2str(avgtime)]);
    fprintf(fid, '\n');
end
fclose(fid);