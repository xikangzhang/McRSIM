% draw shape matrix as 3D blob

figure(1);
membrane(7);
hold on;
quiver3(0,0,0,2,0,0,'b');
quiver3(0,0,0,0,2,0,'k');
quiver3(0,0,0,0,0,2,'r');
hold off;
axis off;