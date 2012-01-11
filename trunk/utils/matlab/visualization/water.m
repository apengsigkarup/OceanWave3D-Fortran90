function c = water(m)
% WATER Linear blue-tone color map.
% WATER(M) return and M-by-3 matrix containing a "water" colormap.
%
%   WATER, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(water)
%
%   Allan P. Engsig-Karup, 15-01-2005. 

if nargin < 1, m = size(get(gcf,'colormap'),1); end
c = min(1,gray(m)*diag([0 1 1-0.5625]));
c(:,3) = c(:,3) + 0.5625;