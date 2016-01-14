function [] = callback_1(str)
f = findobj('tag', str);
s = findobj('tag', 'slider');
a = findobj('parent', f, 'type', 'axes');

temp = get(f, 'userd');
original = temp{1};
frame = original>= get(s, 'value');
set(f, 'userd', {original;frame});
hold on;
imshow(frame, 'parent', a);
hold off;
end
