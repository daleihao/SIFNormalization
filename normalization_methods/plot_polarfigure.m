function plot_polarfigure(VAA,VZA, Refl, cmin, cmax, title_text, istitle, SZA)
n=6;
maxrho = 50;

ContourPolar(VAA,VZA,Refl, n,maxrho, SZA);
%colorbar
colormap(flipud(brewermap(100,'*PiYG')))
caxis([cmin,cmax])
if istitle
title(title_text, 'FontWeight', 'Bold')
end
set(gca, 'linewidth', 1, 'xtick', [], 'fontsize', 12)
end