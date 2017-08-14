% Motion segmentation on Hopkins 12 Real Motion Sequences With Incomplete Data

clear; close all
addpath(genpath('../3rdParty'));
addpath(genpath('../matlab'));
dataPath = '~/research/data/Hopkins155_AdditionalSequences_MissingData';

func = {
    '[missrate, grp, bestRank, minNcutValue] = RSIM_Incomplete(X'', indX'', s, 4, 1);';
    '[missrate, grp, bestRank, minNcutValue] = RSIM_MDD_Incomplete(X'', indX'', s, 4, 1);'
    };

outputFile = 'outputMissing.txt';
fid = fopen(outputFile, 'w');

for funIndex = 1:length(func)
    file = dir(dataPath);
    ii = 0;
    for i = 1:length(file)
        if( (file(i).isdir == 1) && ~strcmp(file(i).name,'.') && ~strcmp(file(i).name,'..') )
            filePath = file(i).name;
            f = dir(fullfile(dataPath, filePath));
            foundValidData = false;
            for j = 1:length(f)
                if( ~isempty(strfind(f(j).name,'_truth.mat')) )
                    ind = j;
                    foundValidData = true;
                    load(fullfile(dataPath, filePath, f(ind).name));
                    break
                end
            end
            
            if(foundValidData)
                N = size(x,2);
                F = size(x,3);
                D = 2*F;
                
                X = reshape(permute(x(1:2,:,:),[1 3 2]),D,N);
                
                % remove short trajectories
                minLen = 5;
                indShort = getShortTrajIndex(m, minLen);
                X(:, indShort) = [];
                m(indShort, :) = [];
                s(indShort) = [];
                indX = kron(m,[1,1])';
                
                X = X.*indX;
                
                rng('default');
                eval(func{funIndex});
                
                ii = ii+1;
                Missrate(ii) = missrate;
                disp([filePath ': ' num2str(100*Missrate(ii)) '%' ', nMotions: ' num2str(max(s)) ', seq: ' num2str(ii)]);
                clear x y s m
            end
        end
    end
    
    avgtol = mean(Missrate);
    medtol = median(Missrate);
    maxtol = max(Missrate);
    stdtol = std(Missrate);
    
    disp('Results on Hopkins 12 Real Motion Sequences With Incomplete Data');
    disp(['Function name: ', func{funIndex}]);
    disp(['Mean: ' num2str(100*avgtol) '%' ', median: ' num2str(100*medtol) '%;'...
        ', max: ' num2str(100*maxtol) '%;' ', std: ' num2str(100*stdtol) '%;']);
    
    fprintf(fid, '%s\n', 'Results on Hopkins 12 Real Motion Sequences With Incomplete Data');
    fprintf(fid, '%s\n', ['Function name: ', func{funIndex}]);
    fprintf(fid, '%s\n', ['Mean: ' num2str(100*avgtol) '%' ', median: ' num2str(100*medtol) '%;'...
        ', max: ' num2str(100*maxtol) '%;' ', std: ' num2str(100*stdtol) '%;']);
    fprintf(fid, '\n');
    
end
fclose(fid);