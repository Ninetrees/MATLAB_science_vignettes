function [ driftStep dUncertainty driftVelocity drift_E dsQuality] = edi_drift_step ( ...
	obsID, ...
	B_t2k, ...
	B_dmpa, ...
	gd_virtual_dmpa, ...
	gd_fv_dmpa, ...
	gd_ID );

	myLibScienceConstants
	global plot_beams rotation_method sinx_wt_Q_xovr use_v10502
% 	disp 'entering edi_drift_step' % debug

	B2n  = norm (B_dmpa, 2);
	Bu   = B_dmpa / B2n;

	switch rotation_method
		% Method 1 ________________________________________________________________
		case 1
			I3    =  [ 1 0 0 ;  0 1 0 ; 0 0 1 ];
			DMPAz = [ 0.0; 0.0; 1.0 ];
			r     = cross (DMPAz, Bu); % rotate DMPAz to B
			rx    = [ ...
			   0.0   r(3) -r(2);
				-r(3)  0.0   r(1);
				 r(2) -r(1)  0.0 ];
			cosTheta = Bu (3); % dot (Bu, DMPAz): ||OCSz|| = ||Bu|| = 1; % Always Bu (3)
			DMPA2BPP  = rx + (cosTheta * I3) + ((1-cosTheta) * (r*r') / sum (r.^2));

		% Method 2 ________________________________________________________________
		case 2
			DMPA2BPPy = [ 0.0; 1.0; 0.0 ];
			DMPA2BPPx = cross (DMPA2BPPy, Bu);
			DMPA2BPPx = DMPA2BPPx / norm (DMPA2BPPx, 2);
			DMPA2BPPy = cross (Bu, DMPA2BPPx);
			DMPA2BPP  = [ DMPA2BPPx'; DMPA2BPPy'; Bu' ]

		% Method 3 ________________________________________________________________
		otherwise
			disp 'rotation method not implemented yet'
			keyboard
	end % switch

	B_bpp = DMPA2BPP * B_dmpa
	% OR
	B_bpp = [ 0.0; 0.0; B2n ];

	% -~-~-~-~-~-~-~-~-~
	% Rotate GDU positions and firing vectors for each GDU into BPP
	gd_virtual_bpp = DMPA2BPP * gd_virtual_dmpa
	gd_fv_bpp      = DMPA2BPP * gd_fv_dmpa

	% -~-~-~-~-~-~-~-~-~
	% Find the most probable beam convergence for the drift step virtual source S*
	% Theory:
	% The beams were rotated into BPP. Those beams are parallel to BPP (perpendicular to B),
	% though shifted in BBPz according to the tilt of the spacecraft in BPP.
	% If we assume a value of zero for the BPPz component of the beam, then we have but to solve
	% for the intersections of the beams in 2D. Here we gather the info that we'll need later:
	% slope and y-intercept. After processing all the beams, we'll calculate the intersections.
	% y = mx + b => m = dy/dx =; b = y - m*x, where x,y = GDU pos

	gd_m_bpp = zeros (1, size (gd_virtual_dmpa, 2), 'double');
	gd_b_bpp = zeros (1, size (gd_virtual_dmpa, 2), 'double');

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
disp (sprintf ('%7.1f ', gd_m_bpp_deg))
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
				gd_virtual_dmpa (:, igd_m_bpp_sorted (iParallel)) = [];
				gd_fv_dmpa (:, igd_m_bpp_sorted (iParallel)) = [];
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

	driftStep     = [NaN; NaN; NaN];
	dUncertainty  = [NaN; NaN; NaN];
	driftVelocity = [NaN; NaN; NaN];
	drift_E       = [NaN; NaN; NaN];

	nBeamIntercepts  = 0;
	dsQuality = 0;
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
% interceptAngleDeg = interceptAngle*rad2deg
% interceptWeights = interceptWeights
		% This test is only valid if NaN is implements above;
		% otherwise, this is always true.
		interceptWeightsSum = nansum (interceptWeights);
		if interceptWeightsSum > 0.0
			% We will check the upper bound, P0=84%, alpha=P1=0.16 ~> +- 0.08 (lower|upper),
			% so div by 2.0, and pass just the upper bound.
			P0 = 0.84;
			edi_stats_alpha = (1.0 - P0) / 2.0;

			% disp 'Two Grubbs tests: x and y'
			ZConfidenceBounds = norminv (1.0 - (1.0 - P0) / 2.0);
			[ GrubbsBeamInterceptMean(1,1), GrubbsBeamInterceptStdDev(1,1), GrubbsBeamIntercepts ] = ...
				edi_drift_step_Grubbs (beamIntercepts (1,:), interceptWeights, edi_stats_alpha);
			ibx = find (isnan (GrubbsBeamIntercepts));

			[ GrubbsBeamInterceptMean(2,1), GrubbsBeamInterceptStdDev(2,1), GrubbsBeamIntercepts ] = ...
				edi_drift_step_Grubbs (beamIntercepts (2,:), interceptWeights, edi_stats_alpha);
			iby = find (isnan (GrubbsBeamIntercepts));
% keyboard
			iIntersectXYOutliers = union (ibx, iby);
			GrubbsBeamIntercepts = beamIntercepts;
			GrubbsBeamIntercepts (:, iIntersectXYOutliers) = [];

			nGrubbsBeamIntercepts          = size (GrubbsBeamIntercepts, 2);
			GrubbsBeamInterceptMean_stdDev = GrubbsBeamInterceptStdDev / sqrt (nGrubbsBeamIntercepts); % x,y mean std dev
% 			disp ( sprintf ('nGrubbsBeamIntercepts %4d : GrubbsInterceptMean, uncertainty: (%+7.3f, %+7.3f) (%+7.3f, %+7.3f)', ...
% 				nGrubbsBeamIntercepts, GrubbsBeamInterceptMean, ZConfidenceBounds*GrubbsBeamInterceptMean_stdDev) ) % debug
% keyboard
			% -~-~-~-~-~-~-~-~-~
			% now we need the drift step...
			virtualSource_bpp = [ GrubbsBeamInterceptMean(1); GrubbsBeamInterceptMean(2); 0.0 ];
			gyroFrequency = (q * B2n * nT2T) / e_mass; % (SI) (|q| is positive here.)
			gyroPeriod    = (twoPi / gyroFrequency);    % (SI) The result is usually on the order of a few ms
% keyboard
			% vE = v in direction of E; T = gyroPeriod
			% ( vE = d/T ) = ExB/|B|^2 ~> d / T * |B|^2 = ExB --- Pacshmann, 1998, 2001, EDI for Cluster
			% Cross from the left with B: B x [] = BxExB
			% where BxExB = E(B dot B) - B(E dot B)
			% The second term is zero because we are assuming E is perpendicular to B.
			% B x [ d/T * |B|^2 = E * |B|^2 ~> E = B x d/T

			% the virtual source S* is the negative of the real drift step; S* is an imaginary point
			% This should be E = B x v, but B, v are swapped here because we need the real drift step (drift velocity),
			% not the virtual source, S*. See relevant publications on Cluster drift step
			% and 'EDI_beams_and_virtual_source_demo_0101.m'.
			E_bpp = cross (virtualSource_bpp, B_bpp) * (1.0e-9 / gyroPeriod); % B_bpp is in nT, and all these calcs are in SI
			driftStep     = -(DMPA2BPP' * virtualSource_bpp);
			dUncertainty  = [ZConfidenceBounds*GrubbsBeamInterceptMean_stdDev; 0.0];
			% Possible future? dUncertainty  = (DMPA2BPP' * [ZConfidenceBounds*GrubbsBeamInterceptMean_stdDev; 0.0])
			driftVelocity = driftStep / gyroPeriod * 1.0e-3; % m/s ~> km/s, per MMS unit standards
			drift_E       = (DMPA2BPP' * E_bpp) * 1.0e3; % convert V/m -> mV/m

			dsQualityWeight = interceptWeightsSum / nGrubbsBeamIntercepts;
			if dsQualityWeight > sinx_wt_Q_xovr(2)
				dsQuality = 3;
			else
				if dsQualityWeight > sinx_wt_Q_xovr(1)
					dsQuality = 2;
				else
					dsQuality = 1;
				end
			end
			% -~-~-~-~-~-~-~-~-~
			if plot_beams
				edi_drift_step_plot ( ...
					obsID, ...
					B_t2k, ...
					gd_ID, ...
					gd_virtual_bpp, ...
					gd_fv_bpp, ...
					B_dmpa, ...
					virtualSource_bpp, ...
					DMPA2BPP, ...
					GrubbsBeamIntercepts, GrubbsBeamInterceptMean, GrubbsBeamInterceptMean_stdDev, ...
					P0, dsQuality);
			end
% 			keyboard
		end % nansum (interceptWeights) > 0.0
	end %  nBeams > 1
end
