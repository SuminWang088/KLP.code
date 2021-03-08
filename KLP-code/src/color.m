 [c,h]=contourf(t1, t2, Z,  'levelstep', 0.08, 'linewidth', 0.5),title('');
map = [0.52 0.8 0.92
0.8 0.7 0.8
0.98 0.94 0.90
1 0.75 0.79
0.94 1 0.94];
colormap(map);
set(h,'ShowText','on');
clabel(c,h,'LabelSpacing',10000);