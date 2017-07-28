function trajR = rotation(traj, theta)

R = [cos(theta), -sin(theta);
    sin(theta), cos(theta)];
trajR = R * traj;

end