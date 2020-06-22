function [ r ] = gkeroSparseR_HDRF( theta )
%  HRDF or albedo vol-kernel (MODIS BRDF/ALBEDO ATBD)
% Angle Input is radian angle
% Code by Dalei Hao,2016/09/21
% Email: dalei.hao.93@gmail.com
%--------------------------------------------------------------------------
  r = -1.284909 -0.166314*theta.^2+0.041840*theta.^3;


