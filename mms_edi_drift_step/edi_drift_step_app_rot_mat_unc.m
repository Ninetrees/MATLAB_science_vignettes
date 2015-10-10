% edi_drift_step_apply_rotation_matrix_uncertainty
% Technically, the rot mat is applied to the col vector as R * v,
% implying that the the first pass of the uncertainty calcs should be
% R (1, 1:3) * v (1:3) ~> v (1:3).
% But MATLAB is faster with col vctors than with row vectors, so the rot_mat_frunc
% passed in should be the transpose.
% There are 2 important !!! points to make:
% 1. Note the inverted rotation and rotation uncertainty matrices.
% 2. Note that the actual rotation is done here. We need the intermediate results,
%    so might as well do it here.

function [rot_vec, rot_vec_unc ] = edi_drift_step_app_rot_mat_unc ( ...
	rot_mat, vec,...
	rot_mat_frunc, vec_frunc);

	R1v_prod_vec = rot_mat (:,1) .* vec;
	R1v_frunc_x  = rot_mat_frunc (:,1) + vec_frunc; % 3D
	R1v_x        = sum (R1v_prod_vec);
	R1v_unc_x    = sum (R1v_frunc_x .* R1v_prod_vec);

	R2v_prod_vec = rot_mat (:,2) .* vec;
	R2v_frunc_y  = rot_mat_frunc (:,2) + vec_frunc; % 3D
	R2v_y        = sum (R2v_prod_vec);
	R2v_unc_y    = sum (R2v_frunc_y .* R2v_prod_vec);

	R3v_prod_vec = rot_mat (:,3) .* vec;
	R3v_frunc_z  = rot_mat_frunc (:,3) + vec_frunc; % 3D
	R3v_z        = sum (R3v_prod_vec);
	R3v_unc_z    = sum (R3v_frunc_z .* R3v_prod_vec);

	rot_vec      = [ R1v_x; R2v_y; R3v_z ];
	rot_vec_unc  = [ R1v_unc_x; R2v_unc_y; R3v_unc_z ];
end