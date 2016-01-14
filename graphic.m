figure('position', [50 50 1000 800]);
axes('unit', 'pix', 'posi', [50 50 696 316]);
imshow(boundary{1});

uicontrol('style', 'slider', 'min', 10000, 'max', 30000, 'value', 22600, ...
    'unit', 'pix', 'position', [800 400 100 20]);
uicontrol('style', 'edit', 'value', 22600, ...
    'unit', 'pix', 'position', [800 380 100 20]);
