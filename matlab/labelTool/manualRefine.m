function loc = manualRefine(vidFrame, loc)

displayPoints(vidFrame, loc)

r = [];
while isempty(r)
    rect = getrect;
    bb = [rect(1), rect(1)+rect(3)-1, rect(2), rect(2)+rect(4)-1];
    isInBbox = (loc(:, 1) >= bb(1) & loc(:, 1) <= bb(2) & loc(:, 2) >= bb(3) & loc(:, 2) <= bb(4));
    loc(isInBbox, :) = [];
    displayPoints(vidFrame, loc);
    r = input('Press Enter to continue. Press q if you are done. \n','s');
end

end

function displayPoints(vidFrame, loc)

imshow(vidFrame);
hold on;
title('Detected features');
plot(loc(:, 1), loc(:, 2), 'g+');
hold off;

end