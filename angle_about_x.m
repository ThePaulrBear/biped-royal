% function angle_out = angle_about_x(angle_in, axis)
% % For a rotation that is purely about the x-axis, this function takes the
% % "angle" and "axis" outputs of the transform sensor and returns the
% % angle rotated about the x-axis. The reason this is needed is that the
% % transform sensor sometimes returns the -x axis so we must switch the sign
% % of the rotation. 
% 
% if angle_in == 0
%     angle_out = 0;
%     return
% end
% 
% assert(abs(axis(1)) == 1, "axis or rotation is not in x direction. axis=(%.2f, %.2f, %.2f)", axis(1), axis(2), axis(3));
% 
% angle_out = angle_in * axis(1);
