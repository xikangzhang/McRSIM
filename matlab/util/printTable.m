function printTable(results, fileName)

if nargin < 2
    fileName = './Table.txt';
end
fid = fopen(fileName,'w');


rowTerm = {'avgtwo', 'medtwo', 'avgthree', 'medthree', 'avgtol', ...
    'medtol', 'time', 'avgtime'};
colTerm = {
        '[missrate,C1, label, gt] = SSCwrapper(X, s);';
    '[missrate, grp, CKSym, index] = SSC_MDD(X, s);';
    '[missrate, grp, index] = SSC_MDD_LA(X, s, camID);';
    '[missrate, label, gt] = LRRwrapper(X, s);';
    '[missrate, label, gt] = LRR_MDD(X, s);';
    '[missrate, grp, index] = LRR_MDD_LA(X, s, camID);';
    '[missrate, grp, bestRank, minNcutValue,W, index] = RSIM(X, s, 2, 1);';
    '[missrate, grp, bestRank, minNcutValue,W, index] = RSIM_MDD(X, s, 2, 1);';
    '[missrate, grp, index] = RSIM_MDD_LA(X, s, 2, 1, camID);';
    '[missrate, grp, bestRank,W, index] = McRSIM(X, s, 2, 1, camID);';
    '[missrate, grp, bestRank,W, index] = McRSIM_MDD(X, s, 2, 1, camID);';
    };

funName = {results.funName};


for i = 1:length(rowTerm)
    for j = 1:length(colTerm)
        ind = strcmp(funName, colTerm{j});
        res = results(ind);
        a = getfield(res, rowTerm{i});
        if j < length(colTerm)
            fprintf(fid, ' %.02f & ', a);
        else
            fprintf(fid, ' %.02f \\\\\n', a);
        end
    end
end

fclose(fid);

end