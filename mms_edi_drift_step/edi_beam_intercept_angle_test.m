% edi_beam_intercept_angle_test.m
% This code should return the smallest angle, 0..90°, twixt two lines that intersect.
% This code returns the abs() of the angle, because that is what edi_drift_step
% needs; otherwise, the sign of the angle reflects the order of subtraction
% in the numerator.
%
% edi_drift_step uses the beam intercept angle to assign a 'confidence' weight
% to the intersection. The maximum confidence is when the beams are perpendicular
% in BPP; the minimum, when the beams are nearly parallel. edi_drift_step does
% not distinguish twixt two beams that intersect, for example, with an
% internal angle of 160° or of 20°. Both receive the same weight.

clc
deg2rad = pi / 180.0;
rad2deg = 180.0 / pi;

for angle = 10: 10: 170
	angle1 = angle * deg2rad;
	angle2 = (180 - angle) * deg2rad;
	slope1 = tan (angle1);
	slope2 = tan (angle2);
	interceptAngle = abs (atan ( (slope2 - slope1) / (1.0 + slope1 * slope2) ) );
	disp (sprintf ('%4d %4d %4.2f %4.2f %6.1f', angle, (180-angle), slope1, slope2, interceptAngle*rad2deg) )
end
