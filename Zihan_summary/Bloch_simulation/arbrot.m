%	function [M] = arbrot(alpha,theta,phi)
%
%	Function returns the rotation matrix M such that
%	completely arbitrary B1 direction with
%   flip angle: alpha in the z-xy plane
%   theta: positive x to positive y
%   phi: z direction
%   all angles in degree
%

% ======================== CVS Log Messages ========================
% by Zihan Ning
% at King's college london 
% 27-May-2025
% ================================================================== 


function [M] = arbrot(alpha,theta,phi)
% ZN: angles in degree

Rx = xrot(alpha);
Ry = yrot(theta);
Rmy = yrot(-theta);
Rz = zrot(phi);
Rmz = zrot(-phi);

M = Rmz*Rmy*Rx*Ry*Rz;
