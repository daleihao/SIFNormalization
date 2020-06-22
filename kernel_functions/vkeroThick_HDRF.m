function [ r ] = vkeroThick_HDRF( theta )
%  HRDF or albedo vol-kernel (MODIS BRDF/ALBEDO ATBD)
% Angle Input is radian angle
% Code by Dalei Hao,2016/09/21
% Email: dalei.hao.93@gmail.com
%--------------------------------------------------------------------------

  r = -0.007574-0.070987*theta.^2+0.307588*theta.^3;




