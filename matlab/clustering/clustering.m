function clustering(X, n_classes, gt, actions, opt)

% gt = action_labels(action_ind)';
% [label,X_centers,sD] = ncutJLD(X,n_classes,opt);
[label,X_centers,sD] = kmeansJLD(X,n_classes,opt);
v = perms(actions);
acc = zeros(1,size(v,1));
for i = 1:length(acc)
    acc(i) = nnz(v(i,label)==gt)/length(gt);
end
[accuracy,ind] = max(acc);
accuracy
label = v(ind,label);
confusion_matrix = zeros(n_classes, n_classes);
for i = 1:n_classes
    temp = find(gt == actions(i));
    confusion_matrix(i, :) = hist(label(temp), actions) / length(temp);
end

end