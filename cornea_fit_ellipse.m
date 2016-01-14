function [ pt_ellipse ] = cornea_fit_ellipse( cornea )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

A = EllipseDirectFit( cornea );
a = A(1); b = A(2); c = A(3); d = A(4); e = A(5); f = A(6);

ymax = max(cornea(:,2)); ymin = min(cornea(:,2));

y = (ymin:ymax)';

square_root = sqrt( (b*y+d).^2 - 4*a*(c*y.^2 + e*y+f) );

x = (0.5/a) * (-b*y-d - square_root );

pt_ellipse = [x,y];


end

