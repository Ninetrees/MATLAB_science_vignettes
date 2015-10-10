function [ ...
		d_bpp, ...
		d_SD_bpp, ...
		d_CI_bpp, ...
		d_sdcs, ...
		d_SD_sdcs, ...
		d_CI_sdcs, ...
		d_quality, ...
		v_sdcs, ...
		E_sdcs, ...
		E_sdcs_unc ...
	] = edi_drift_step ( ...
		obsID, ...
		B_t2k, ...
		B_sdcs, ...
		B_SD_sdcs, ...
		gd_virtual_sdcs, ...
		gd_fv_sdcs, ...
		gd_ID ...
	);

	myLibScienceConstants
	global plot_beams rot_mat_unc rotation_method sinx_wt_Q_xovr use_v10502
% 	disp 'entering edi_drift_step' % debug

	B2n  = norm (B_sdcs, 2);
	% The first step in propagating uncertainty to the E-field starts with B2n.
	% The 2norm is the sqrt of the sum of the squares.
	% Rule 1: Note use of fractional uncertainties.
	% Rule 2: Uncertainties for sums: sum uncertainties (_unc) ordinary or in ***quadrature.
	% Rule 3: Uncertainties for products: sum fractional uncertainties (_frunc).
	% Rule 4: If q = x^n, then q_frunc = |n| * x_frunc.
	% Rule 5: If q = Cx, where C has no uncertainty, q_unc = |C| * x_unc (NOT _frunc).
	% Rule 6: Unc in functions of one variable: q_unc = |dq/dx| * x_unc.
	% Rule 7: Use SDOM where available. [1] Sect. 4.4.
	% Rule 8: Convert _frunc to _unc: q_unc = |q| * q_frunc
	% ***Because the goal is to use fewer beams per set, no quadrature is used.

	% Calc uncertainty for B2n. These lines are all about the UNCERTAINTY, NOT B2n.
	% Here each Bx|y|z is squared, then those squares are summed, then sqrt (sum).
	% The uncertainty: Rule 1, 4, 2, 4.
	B_frunc = abs (B_SD_sdcs ./ B_sdcs); % Rule 1
	% Techincally we need || here for B?_unc, but we are squaring them below, so could ignore.
	Bx_unc  = abs ((2.0 * B_frunc (1)) * B_sdcs (1)); % Rule 4. Bx|y|z is raised to a power.
	By_unc  = abs (2.0 * B_frunc (2) * B_sdcs (2)); % We need unc, not frunc, for Rule 2
	Bz_unc  = abs (2.0 * B_frunc (3) * B_sdcs (3));
	% It is tempting to combine the next two into 1 line, but we need B2n_frunc later.
	B2n_frunc = 0.5 * ((Bx_unc + By_unc + Bz_unc) / B2n); % Rule 2, 1, 4
	B2n_unc = B2n_frunc * B2n; % Total uncertainty in B2n
% 	disp 'B stats' % V&V
% 	[ B_sdcs B_SD_sdcs [B2n_unc;0;0] ] % V&V

	Bu_sdcs       = B_sdcs / B2n; % Bu has uncertainty
	Bu_sdcs_frunc = B_frunc + B2n_frunc; % 3D
	Bu_sdcs_unc   = abs (Bu_sdcs_frunc .* Bu_sdcs); % 3D

	% We have not decided how to handle the uncertainty of the rotation matrix, 2015-10-01
	% So these uncertainty calcs are wrapped in rot_mat_unc flag checks
	switch rotation_method
		% Method 1 ________________________________________________________________
		case 1
			I3    =  [ 1 0 0 ;  0 1 0 ; 0 0 1 ];
			DMPAz = [ 0.0; 0.0; 1.0 ];
			r     = cross (DMPAz, Bu_sdcs); % rotate DMPAz to B
			rx    = [ ...
			   0.0   r(3) -r(2);
				-r(3)  0.0   r(1);
				 r(2) -r(1)  0.0 ];
			cosTheta = Bu_sdcs (3); % dot (Bu, DMPAz): ||OCSz|| = ||Bu|| = 1; % Always Bu (3)
			SDCS2BPP = rx + (cosTheta * I3) + ((1-cosTheta) * (r*r') / sum (r.^2))

		% Method 2 ________________________________________________________________
		case 2
			SDCS2BPPy   = [ 0.0; 1.0; 0.0 ];
			% SDCS2BPPx = cross (SDCS2BPPy, Bu_sdcs); % Replace with next line shortcut
			SDCS2BPPx   = [ Bu_sdcs(3); 0.0; -Bu_sdcs(1) ]; % No new uncertainty to here, just that of Bu
			S2Bx2n      = norm (SDCS2BPPx, 2);
			% -~-~-~-~-~-~-~-~-~
			if rot_mat_unc % uncertainty for S2Bx2n
				S2Bx_frunc   = abs ([ Bu_sdcs_frunc(3); 0.0; Bu_sdcs_frunc(1) ]);
				S2Bxx_unc    = abs (2.0 * Bu_sdcs_frunc (3) * Bu_sdcs (3));
				S2Bxy_unc    = 0.0;
				S2Bxz_unc    = abs (2.0 * Bu_sdcs_frunc (1) * Bu_sdcs (1));
				S2Bx2n_frunc = 0.5 * (S2Bxx_unc + S2Bxz_unc) / S2Bx2n; %  skip S2Bxy_unc = 0.0
				S2Bx2n_unc   = S2Bx2n_frunc * S2Bx2n; % always non-negative
			end
			% -~-~-~-~-~-~-~-~-~
			SDCS2BPPx = SDCS2BPPx / S2Bx2n; % Normalize SDCS2BPPx to unit vector
			% -~-~-~-~-~-~-~-~-~
			if rot_mat_unc % uncertainty for SDCS2BPPx normalization
				% Bu_sdcs_frunc is a surrogate for the unnormalized SDCS2BPPx
				% S2Bx_frunc = [ Bu_sdcs_frunc(3); 0.0; Bu_sdcs_frunc(1) ];
				S2Bx_frunc = S2Bx_frunc + S2Bx2n_frunc;
				S2Bx_unc   = S2Bx_frunc .* SDCS2BPPx; % 3D
			end
			% -~-~-~-~-~-~-~-~-~
			SDCS2BPPy = cross (Bu_sdcs, SDCS2BPPx);
			% -~-~-~-~-~-~-~-~-~
			if rot_mat_unc % uncertainty for the cross product
				% Each dim is the same set of calcs, but with different uncertainties: 2 products and 1 sum
				% Note that while the cross product has a sign change in j, the fruncs do not
				S2By_ifrunc1 = Bu_sdcs_frunc(2) + S2Bx_frunc(3);
				S2By_ifrunc2 = Bu_sdcs_frunc(3) + S2Bx_frunc(2);
				S2By_jfrunc1 = Bu_sdcs_frunc(1) + S2Bx_frunc(3);
				S2By_jfrunc2 = Bu_sdcs_frunc(3) + S2Bx_frunc(1);
				S2By_kfrunc1 = Bu_sdcs_frunc(1) + S2Bx_frunc(2);
				S2By_kfrunc2 = Bu_sdcs_frunc(2) + S2Bx_frunc(1);
				% Calculate _unc from _frunc, quad add
				S2By_unc (1) = S2By_ifrunc1 * SDCS2BPPy(1) + S2By_ifrunc2 * SDCS2BPPy(1);
				S2By_unc (2) = S2By_jfrunc1 * SDCS2BPPy(2) + S2By_jfrunc2 * SDCS2BPPy(2);
				S2By_unc (3) = S2By_kfrunc1 * SDCS2BPPy(3) + S2By_kfrunc2 * SDCS2BPPy(3);
			end
			% -~-~-~-~-~-~-~-~-~
			SDCS2BPP   = [ SDCS2BPPx'; SDCS2BPPy'; Bu_sdcs' ];
			S2By_frunc = S2By_unc' ./ SDCS2BPPy;
			SDCS2BPP_frunc = [ ...
				[ S2Bx_frunc(1) S2Bx_frunc(2) S2Bx_frunc(3) ]
				[ S2By_frunc(1) S2By_frunc(2) S2By_frunc(3) ]
				[ Bu_sdcs_frunc(1)   Bu_sdcs_frunc(2)   Bu_sdcs_frunc(3)   ] ]; % V&V
			SDCS2BPP_unc = [ ...
				[ S2Bx_unc(1)    S2Bx_unc(2)    S2Bx_unc(3) ]
				[ S2By_unc(1)    S2By_unc(2)    S2By_unc(3) ]
				[ Bu_sdcs_unc(1) Bu_sdcs_unc(2) Bu_sdcs_unc(3) ] ];

		% Method 3 ________________________________________________________________
		otherwise
			disp 'Specified rotation method not implemented yet.'
			keyboard
	end % switch
	BPP2SDCS       = SDCS2BPP';
	BPP2SDCS_unc   = SDCS2BPP_unc';
	BPP2SDCS_frunc = SDCS2BPP_frunc';

	% B_bpp = SDCS2BPP * B_sdcs; % V&V
	% OR
	B_bpp = [ 0.0; 0.0; B2n ];

	% -~-~-~-~-~-~-~-~-~
	% Rotate GDU positions and firing vectors for each GDU into BPP
	gd_virtual_bpp = SDCS2BPP * gd_virtual_sdcs;
	gd_fv_bpp      = SDCS2BPP * gd_fv_sdcs;

	% -~-~-~-~-~-~-~-~-~
	% Find the most probable beam convergence for the drift step virtual source S*
	% Theory:
	% The beams were rotated into BPP. Those beams are parallel to BPP (perpendicular to B),
	% though shifted in BBPz according to the tilt of the spacecraft in BPP.
	% If we assume a value of zero for the BPPz component of the beam, then we have but to solve
	% for the intersections of the beams in 2D. Here we gather the info that we'll need later:
	% slope and y-intercept. After processing all the beams, we'll calculate the intersections.
	% y = mx + b => m = dy/dx =; b = y - m*x, where x,y = GDU pos

	gd_m_bpp = zeros (1, size (gd_virtual_sdcs, 2), 'double');
	gd_b_bpp = zeros (1, size (gd_virtual_sdcs, 2), 'double');

	% -~-~-~-~-~-~-~-~-~
	% A gd12 beam originates in gun 1 and is detected in det 2.
	% Find all the beam slopes and y-intercepts
	gd_m_bpp = gd_fv_bpp (2,:) ./ gd_fv_bpp (1,:);
	gd_b_bpp = gd_virtual_bpp (2,:) - gd_m_bpp.* gd_virtual_bpp (1,:);

	% v1.05.xx, new calc: find beams that are parallel within 4° (slope ~0.07)
	% This REMOVES beams that are nearly parallel
	if use_v10502
		iParallel = [0];
		while ~isempty (iParallel)
			gd_m_bpp_deg = atand (gd_m_bpp); % slopes in degrees

			% sort the slopes from lowest to highest, mono
			[gd_m_bpp_sorted, igd_m_bpp_sorted] = sort (gd_m_bpp_deg);
			[gd_m_bpp_sorted; gd_m_bpp_deg(igd_m_bpp_sorted)];
			% Calc the diff twixt adjacent slopes... not the best, but a prototype.
			% Some of the deltas may be negative; make them all non-negative.
			gd_m_bpp_diff = abs (diff (gd_m_bpp_sorted));
			% This is arbitrary: at what intersection angle do we want to penalize beams?
			% If the max beam uncertainty is ~+-2°, and 2 beams are fired parallel,
			% then the largest angle twixt them, due to uncertainty, is 4°.
			% IOW, if two beams are within 4° of each other, we are not certain that
			% they are not supposed to be parallel.
			iParallel = find (gd_m_bpp_diff < 4);

			% This is a prototype. It simply tosses all but the first beam in a series
			% of beams that meet the parallelism test above. But the side effect is that
			% the third beam in the series, which might be /more than/ 4° from the first,
			% gets deleted if it is within 4° of the second.
			if ~isempty (iParallel)
				disp ( sprintf ('n iParallel: %2d of %2d\n', length(iParallel), length(gd_m_bpp)) )
				gd_m_bpp (igd_m_bpp_sorted (iParallel)) = [];
				gd_b_bpp (igd_m_bpp_sorted (iParallel)) = [];
				gd_virtual_sdcs (:, igd_m_bpp_sorted (iParallel)) = [];
				gd_fv_sdcs (:, igd_m_bpp_sorted (iParallel)) = [];
				gd_ID (igd_m_bpp_sorted (iParallel)) = [];
				gd_virtual_bpp (:, igd_m_bpp_sorted (iParallel)) = [];
				gd_fv_bpp (:, igd_m_bpp_sorted (iParallel)) = [];
			end
		end
	end
% keyboard
	% Find the S* in BPP, using BPP FV convergence
	% preAlloc beamIntercepts based on nBeams: (nBeams - 1) * nBeams / 2
	% -~-~-~-~-~-~-~-~-~
	nBeams           = uint32 (length (gd_m_bpp));
	nBeamIntercepts  = (nBeams-1) * nBeams / 2;
	beamIntercepts   = zeros (2, nBeamIntercepts, 'double');
	interceptWeights = zeros (1, nBeamIntercepts, 'double');

	d_bpp    = [NaN; NaN; NaN];
	d_CI_bpp = [NaN; NaN; NaN];
	v_sdcs   = [NaN; NaN; NaN];
	E_sdcs   = [NaN; NaN; NaN];

	nBeamIntercepts  = 0;
	d_quality = 0;
	if nBeams > 1 % can't have an intersection with just 1 beam

		% find the intercepts for all beams
		for i = 1: nBeams-1
			for j = i+1: nBeams
				XY = [ ...
					-gd_m_bpp(i) 1 ;
					-gd_m_bpp(j) 1 ];
				b = [ gd_b_bpp(i) ; gd_b_bpp(j) ];
				nBeamIntercepts = nBeamIntercepts + 1;
				beamIntercepts (:, nBeamIntercepts) = XY \ b;
			end
		end

		% Not determined as of v10501 is what angles really get what weight
		% We tested sin^2, and got some satisfaction, but we really want to let most beams
		% count, and sharply penalize beams that intersect at small angles
		% changed to sin^4
		% The check for macroBeamCheckAngle removes intercepts from CALCS if the
		% intercept angle is less than macroBeamCheckAngle.

		macroBeamCheckAngle = atan(tand(5));
		nBeamIntercepts = 0;
		for i = 1: nBeams-1
			for j = i+1: nBeams
				nBeamIntercepts = nBeamIntercepts + 1;

				interceptAngle (nBeamIntercepts) = abs (atan ( ...
					(gd_m_bpp(j) - gd_m_bpp(i)) / ...
					(1.0 + gd_m_bpp(i) * gd_m_bpp(j)) ) );

% 				if (abs (interceptAngle (nBeamIntercepts)) > macroBeamCheckAngle) % v1.05.01
					interceptWeights (1, nBeamIntercepts) = sin (interceptAngle (nBeamIntercepts))^4;
% 				else
% 					interceptWeights (1, nBeamIntercepts) = NaN;
% 				end
			end
		end

		% This test is only valid if NaN is implemented above;
		% otherwise, this is always true.
		interceptWeightsSum = nansum (interceptWeights);
		if interceptWeightsSum > 0.0
			% We will check the upper bound, P0=84%, alpha=P1=0.16 ~> +- 0.08 (lower|upper),
			% so div by 2.0, and pass just the upper bound.
			% See 'Accuracy, Error, Precision, and Uncertainty.txt'
			P0 = 0.84;
			edi_stats_alpha = (1.0 - P0)  / 2.0;
			z_score = norminv (1 - edi_stats_alpha);
			% disp 'Two Grubbs tests: x and y'
			[ GrubbsBeamInterceptMean(1,1), GrubbsBeamInterceptSD(1,1), GrubbsBeamIntercepts ] = ...
				edi_drift_step_Grubbs (beamIntercepts (1,:), interceptWeights, edi_stats_alpha);
			ibx = find (isnan (GrubbsBeamIntercepts));

			[ GrubbsBeamInterceptMean(2,1), GrubbsBeamInterceptSD(2,1), GrubbsBeamIntercepts ] = ...
				edi_drift_step_Grubbs (beamIntercepts (2,:), interceptWeights, edi_stats_alpha);
			iby = find (isnan (GrubbsBeamIntercepts));
% keyboard
			iIntersectXYOutliers = union (ibx, iby);
			GrubbsBeamIntercepts = beamIntercepts;
			GrubbsBeamIntercepts (:, iIntersectXYOutliers) = [];

			nGrubbsBeamIntercepts   = size (GrubbsBeamIntercepts, 2);
			% The uncertainty in S_star_bpp, /// at the 68% conficence level ///
			% Later, we will adjust this value to a different confidence level
			GrubbsBeamInterceptSDOM = GrubbsBeamInterceptSD / sqrt (nGrubbsBeamIntercepts);

% 			disp ( sprintf ('nGrubbsBeamIntercepts %4d : GrubbsInterceptMean, uncertainty: (%+7.3f, %+7.3f) (%+7.3f, %+7.3f)', ...
% 				nGrubbsBeamIntercepts, GrubbsBeamInterceptMean, z_score*GrubbsBeamInterceptSDOM) ) % debug
% keyboard

			% -~-~-~-~-~-~-~-~-~
			% now we need the drift step...
			% gyroFrequency = (q * B2n * nT2T) / mass_e; % (SI) |q| is positive here.
			gyroFrequency     = q_over_mass_e_nT2T * B2n;     % (SI) |q| is positive here.
			gyroFrequency_unc = q_over_mass_e_nT2T * B2n_unc; % Rule 5.
			gyroPeriod        = twoPi / gyroFrequency;        % (SI) Usually on the order of a few ms
			gyroPeriod_frunc  = twoPi * gyroFrequency_unc / gyroFrequency; % Rule1, 5.
			gyroPeriod_unc    = gyroPeriod_frunc * gyroPeriod;
% 			disp 'gyroFreq gyroPeriod stats'
% 			[ gyroFrequency gyroFrequency_unc gyroPeriod gyroPeriod_unc ] % V&V

			% vE = v in direction of E; T = gyroPeriod
			% ( vE = d/T ) = ExB/|B|^2 ~> d / T * |B|^2 = ExB --- Pacshmann, 1998, 2001, EDI for Cluster
			% Cross from the left with B: B x [] = BxExB
			% where BxExB = E(B dot B) - B(E dot B)
			% The second term is zero because we are assuming E is perpendicular to B.
			% B x [ d/T * |B|^2 = E * |B|^2 ~> E = B x d/T

			% the virtual source S* is the negative of the real drift step; S* is an imaginary point
			S_star_bpp  = [ GrubbsBeamInterceptMean; 0.0 ];
			d_bpp       = -S_star_bpp; % Note the minus sign
			d_SD_bpp    = [GrubbsBeamInterceptSD; 0.0];   % > non-neg by def
			d_SDOM_bpp  = [GrubbsBeamInterceptSDOM; 0.0]; % ditto
			d_bpp_frunc = d_SDOM_bpp(1:2) ./ d_bpp(1:2); % Rule 1, 7.
			d_bpp_frunc (3) = 0.0;
			% Compute uncertainty (confidence interval CI) = z_score * SDOM
			d_CI_bpp = z_score * d_SDOM_bpp;

			% d_bpp is assumed to be zero in BPPz, so the uncertainty is just in x,y
% 			disp 'd stats'
% 			[ d_bpp d_SD_bpp d_SDOM_bpp d_CI_bpp ] % V&V

			gyroPeriod_1e3 = gyroPeriod * 1.0e3;
			v_bpp = d_bpp / gyroPeriod_1e3; % m/s ~> km/s, per MMS unit standards
			v_bpp_frunc = (d_bpp_frunc(:) + gyroPeriod_frunc) * 1.0e-3; % Rule 3, 5.
% 			v_bpp_unc = v_bpp_frunc .* v_bpp; % V&V
% 			disp 'v stats'
% 			[ v_bpp v_bpp_frunc ] % V&V

% 			strB_time = datestr (spdftt2000todatenum (B_t2k), 'HH:MM:ss'); % V&V
% 			disp (sprintf('84%% d_BPP confidence intervals %9s (x, y)= ( %+6.2f, %+6.2f )', ...
% 				strB_time, d_CI_bpp(1), d_CI_bpp(2))) % V&V

			% B_bpp is in nT, and all these calcs are in SI
			% v_sdcs above is in /// km/s ///
			% convert nT -> T, km/s -> m/s, V/m -> mV/m : 1.0e-9 * 1.0e3 * 1.0e3 ~> 1.0e-3
			% This should be E = B x v = B x d/t, not the virtual source S*.
			% See relevant publications on Cluster drift step
			% and 'EDI_beams_and_virtual_source_demo_0101.m'.
			% B_bpp = (0, 0, Bz=B2n); v_bpp = (vx, vy, 0)
			% cross (B_bpp, v_bpp) ~> B2n * [ -v_bpp(2); v_bpp(1); 0.0 ]
			E_bpp       = [ -v_bpp(2); v_bpp(1); 0.0 ] * (B2n * 1.0e-3);
			E_bpp_frunc = [ v_bpp_frunc(2)+B2n_frunc; v_bpp_frunc(1)+B2n_frunc; 0.0 ] * 1.0e-3; % Rule 3, 5.
			E_bpp_unc   = E_bpp_frunc .* E_bpp;

% 			disp 'E stats' % V&V
% 			[ E_bpp E_bpp_frunc E_bpp_unc ] % V&V

			% !!! See important help note in the edi_drift_step_app_rot_mat_unc header.
			% It explains why SDCS2BPP replaces BPP2SDCS.
			% d_sdcs = BPP2SDCS * d_bpp;
			[ d_sdcs, d_sdcs_unc ] = edi_drift_step_app_rot_mat_unc (SDCS2BPP, d_bpp, SDCS2BPP_frunc, d_bpp_frunc);

			d_SD_sdcs  = abs (BPP2SDCS * d_SD_bpp);
			d_CI_sdcs  = abs (BPP2SDCS * d_CI_bpp);

			v_sdcs = d_sdcs / gyroPeriod_1e3; % m/s ~> km/s, per MMS unit standards
			v_sdcs_frunc = (d_sdcs_unc ./ d_sdcs) + gyroPeriod_frunc;
			v_sdcs_unc = v_sdcs_frunc .* v_sdcs;

			% E_sdcs     = BPP2SDCS * E_bpp;
			% E_sdcs_unc = abs (BPP2SDCS * E_bpp_unc)
			[ E_sdcs, E_sdcs_unc ] = edi_drift_step_app_rot_mat_unc (SDCS2BPP, E_bpp, SDCS2BPP_frunc, E_bpp_frunc);

% keyboard

			d_qualityWeight = interceptWeightsSum / nGrubbsBeamIntercepts;
			if d_qualityWeight > sinx_wt_Q_xovr(2)
				d_quality = 3;
			else
				if d_qualityWeight > sinx_wt_Q_xovr(1)
					d_quality = 2;
				else
					d_quality = 1;
				end
			end
			% -~-~-~-~-~-~-~-~-~
			if plot_beams
				edi_drift_step_plot ( ...
					obsID, ...
					B_t2k, ...
					B_sdcs, ...
					gd_ID, ...
					gd_virtual_bpp, ...
					gd_fv_bpp, ...
					SDCS2BPP, ...
					S_star_bpp, d_SDOM_bpp, d_CI_bpp, ... % d_SDOM_bpp == S_star_SDOM_bpp
					GrubbsBeamIntercepts, ...   % for plotting dot on counted intersections
					P0, z_score, d_quality);
			end
% 			keyboard
		end % nansum (interceptWeights) > 0.0
	end %  nBeams > 1
end
