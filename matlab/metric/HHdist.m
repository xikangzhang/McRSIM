function D = HHdist(HH1,HH2,opt)

isSymmetric = false;
if isempty(HH2)
    HH2 = HH1;
    isSymmetric = true;
end
if ~isSymmetric
    m = length(HH1);
    n = length(HH2);
    D = zeros(m,n);
    for i = 1:m
        for j = 1:n
            if strcmp(opt.metric,'JBLD')
                D(i,j) = JBLD(HH1{i},HH2{j});
            elseif strcmp(opt.metric,'JBLD_denoise')
                D(i,j) = gramDist_cccp(pcaClean(HH1{i}), HH2{j}, opt);
            elseif strcmp(opt.metric,'binlong')
                D(i,j) = 2 - norm(HH1{i}+HH2{j},'fro');
            elseif strcmp(opt.metric,'AIRM')
                D(i,j) = AIRM(HH1{i},HH2{j});
            elseif strcmp(opt.metric,'LERM')
                D(i,j) = LERM(HH1{i},HH2{j});
            elseif strcmp(opt.metric,'KLDM')
                D(i,j) = KLDM(HH1{i},HH2{j});
            elseif strcmp(opt.metric,'SubspaceAngleFast')
                D(i,j) = subspace(HH1{i},HH2{j});
            elseif strcmp(opt.metric,'SubspaceAngle')
                D(i,j) = SubspaceAngle(HH1{i},HH2{j},opt.SA_thr);
            end
        end
    end
else
    m = length(HH1);
    D = zeros(m);
    for i = 1:m
        for j = i:m
            if strcmp(opt.metric,'JBLD')
                D(i,j) = JBLD(HH1{i},HH2{j});
            elseif strcmp(opt.metric,'JBLD_denoise')
                D(i,j) = gramDist_cccp(pcaClean(HH1{i}), HH2{j}, opt);
            elseif strcmp(opt.metric,'binlong')
                D(i,j) = 2 - norm(HH1{i}+HH2{j},'fro');
            elseif strcmp(opt.metric,'AIRM')
                D(i,j) = AIRM(HH1{i},HH2{j});
            elseif strcmp(opt.metric,'LERM')
                D(i,j) = LERM(HH1{i},HH2{j});
            elseif strcmp(opt.metric,'KLDM')
                D(i,j) = KLDM(HH1{i},HH2{j});
            elseif strcmp(opt.metric,'SubspaceAngleFast')
                D(i,j) = subspace(HH1{i},HH2{j});
            elseif strcmp(opt.metric,'SubspaceAngle')
                D(i,j) = SubspaceAngle(HH1{i},HH2{j},opt.SA_thr);
            end
        end
    end
    D = D + D';
end
    

end