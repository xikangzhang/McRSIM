function displayTraj(y, index, verbose)

if nargin < 3
    figure;
end

hold on;
for j=1:length(index)
    t = squeeze(y(1:2,index(j),:));
    plot(t(1,:), t(2,:), '*-');
end
hold off;

end