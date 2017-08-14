% Motion segmentation on Hopkins 155 with gross contamination

clear; close all
addpath(genpath('../3rdParty'));
addpath(genpath('../matlab'));
dataPath = '~/research/data/Hopkins155';

func = {
    '[missrate, grp, bestRank, minNcutValue,W, index] = RSIM(X, s);';
    '[missrate, grp, bestRank, minNcutValue,W, index] = RSIM_MDD(X, s, 4, 1, 1);';
    };

outputFile = 'outputGross.txt';
fid = fopen(outputFile, 'w');

for funIndex = 1:length(func)
    rng('default');
    file = listFolder(dataPath);
    ii = 0;
    ii2 = 0;
    ii3 = 0;
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
                A = rand(D, N);
                X(indX==0) = A(indX==0);
                
                eval(func{funIndex});
                
                if ~exist('label','var') && ~exist('gt','var')
                    nCluster = size(grp, 2);
                    label = grp * (1:nCluster)';
                    gt = index(s)';
                end
                tc = nnz(label~=gt);
                clear label gt
                
                ii = ii+1;
                Missrate(ii) = missrate;
                TC(ii) = tc;
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
    end
    time = toc;
    avgtime = time/ii;
    
    avgtol = mean(Missrate);
    medtol = median(Missrate);
    avgtwo = mean(Missrate2);
    medtwo = median(Missrate2);
    avgthree = mean(Missrate3);
    medthree = median(Missrate3);
    
    disp('Results on Hopkins155 with gross corrupted entries')
    disp(['Mean of two: ' num2str(100*avgtwo) '%' ', median of two: ' num2str(100*medtwo) '%;']);
    disp(['Mean of three: ' num2str(100*avgthree) '%' ', median of three: ' num2str(100*medthree) '%;']);
    disp(['Mean of all: ' num2str(100*avgtol) '%' ', median of all: ' num2str(100*medtol) '%;']);
    
    fprintf(fid, '%s\n', 'Results on Hopkins155, half rotated 45 degree, with gross corrupted entries');
    fprintf(fid, '%s\n', ['Function name: ', func{funIndex}]);
    fprintf(fid, '%s\n', ['Mean of two: ' num2str(100*avgtwo) '%' ', median of two: ' num2str(100*medtwo) '%;']);
    fprintf(fid, '%s\n', ['Mean of three: ' num2str(100*avgthree) '%' ', median of three: ' num2str(100*medthree) '%;']);
    fprintf(fid, '%s\n', ['Mean of all: ' num2str(100*avgtol) '%' ', median of all: ' num2str(100*medtol) '%;']);
    fprintf(fid, '%s\n', ['total running time: ', num2str(time), ', average running time per sequence: ', num2str(avgtime)]);
    fprintf(fid, '\n');
    
end
fclose(fid);