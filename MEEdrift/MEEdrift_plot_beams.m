% MEEdrift_plot_beams
% disp 'MEEdrift_plot_beams ...' % V&V

%{
	We come here with iiSorted_beam_t2k set
	1. Draw the spacecraft plane rotated into BPP
	2. Draw 3D axis for reference
	3. For each beam data point
		 a. Draw the center GDU
		 b. Draw GDU locations and beams for 4 beams either side (in ssm) of this one, AND within 2 s of this beam
		 c. Find the 2 EFW data points whose times (in ssm) bracket the center beam time (ssm)
		    Interpolate E-field values twixt them wrt time
		 c. Plot the virtual source S* for each beam as a vector from the origin
		    (the detector of the virtual spacecraft) to S*

1. Rotates a chosen EDI beam into BPP, along with +- 4..8 beams on either side, based on a time window.
2. Ditto the GDUs.
3. Calculates the target based on beam convergence.
4. Interpolates 2D SDP data to the center EDI beam. Uses E.B=0 to create ((E_EdotBeq0)xB)_driftStep.
5. Plots the GDUs, EDI beams, unit vectors for B, E_EdotBeq0, and the driftStep, in, say, BPP.
   The point of the unit vectors to show that it all makes sense.
%}
format short g  % +, bank, hex, long, long g, rat, short, short g, short eng
fclose all;     % nothing should be open, but if an error occurs during reentrancy, this should handle any hanging open file handles

beamIntercepts     = zeros (2, (9*10)/2, 'double'); % each beam intersects every other beam; assumes center beam + 4 on either side
edp_E_EdB_BPP      = zeros (3, 1, 'double'); % pre-allocate for debugging
edp_E_interp_d_dmpa= zeros (3, 1, 'double');
GDU_Loc            = zeros (3, 2, nEDP, 'double'); % In BPP, X & Y, Z=0; 1 vector for each GDU
nbeamIntercepts    = uint32 (0);
nTargetBeams       = uint32 (0);
targetBeamsGDU_Loc = zeros (3, 20, 'double'); % needs to be high enough for max expected beam count
targetBeamAngle    = zeros (3, 20, 'double'); % needs to be high enough for max expected beam count
targetBeam2DSlope  = zeros (1, 20, 'double');
targetBeam2Db      = zeros (1, 20, 'double');

if Use_OCS_plot
	set (0, 'CurrentFigure', hDMPA_plot) % hDMPA_plotElements
	clf
else
	set (hDMPA_plot, 'visible', 'off');
end

set (0, 'CurrentFigure', hBPP_plot) % hBPP_plotElements DMPA2BPP
clf
% set (hBPP_plot, 'visible', 'off');
disp '---'

% iCenterBeam is an index, not a time; points to the center beam used for calculations
iCenterBeam = iSorted_beam_t2k (iiSorted_beam_t2k);
% index to B_dmpa B-field record that matches this CenterBeam EDI record
edi_BdvE_recnum = edi_xref2_BdvE (iCenterBeam);
iedi_BdvE_dmpa  = find (BdvE_xref2_edi == edi_BdvE_recnum);

centerBeamB    = edi_B_dmpa (:, iedi_BdvE_dmpa);
centerBeamB2n  = norm (centerBeamB, 2);
centerBeamB_u  = centerBeamB / centerBeamB2n; % unit vector
centerBeamE    = edi_E_dmpa (:, iedi_BdvE_dmpa);
centerBeamEu   = centerBeamE / norm (centerBeamE, 2);
centerBeam_dn  = edi_gd_beam_dn (iCenterBeam);
centerBeam_t2k = edi_gd_beam_t2k (iCenterBeam);
disp ([ 'Center beam time:  ', char(spdfencodett2000(centerBeam_t2k)) ]) % V&V

% virtual source S* == -(drift step) :: drift step ~> drift velocity ~> drift E
% plot S*, but use drift step, drift E for calcs
centerBeam_d  = edi_d_dmpa (:, iedi_BdvE_dmpa);
centerBeam_du = centerBeam_d / norm (centerBeam_d, 2);
S_star_dmpa   = -centerBeam_d;
S_star_dmpa_u = -centerBeam_du;

% This note must remain in all versions of MEEdrift. There are two important
% concepts here regarding coordinate systems.
%
% BCS, for example, is a natural
% surrogate for GSE. In a 3D view of BCS, we see a complete view of the [expanded]
% spacecraft, the GDUs in their positions as they detect coded electron beams,
% the 3D firing angles of the beams, E-, B-, and ExB-normal vectors representing those fields,
% and indication of the drift step. It is natural to view this plot from a point
% above the xy plane. A standard 3D plot shows +x pointing toward the viewer's left,
% +y pointing somewhat to the viewer's right, and +z pointing up.
%
%   When we want to view the plot in BPP, we imagine that we are still in BCS,
% looking "down" on the spacecraft. The problem comes when we do the math to rotate
% elements from DMPA to BPP. If B is pointing in the +DSIz direction, all is well.
% If B is pointing in the -DSIz direction, the spacecraft outline, GDUs, and beams
% all get rotated "upside down". There is no loss of generality, but mentally
% "swapping" views as B changes from -DSIz to +DSIz can be uncomfortable. So when
% B points toward -DSIz, we "invert" it for the express purpose of plotting a
% consistent view of the spacecraft outline, the GDUs, and the beams. The rest of
% the plotted values remain as they were, and --- MORE IMPORTANT --- the math does
% not change. We still calculate the real drift step, based on real data in DMPA.
%
%   As a result of inverting B and calculating a rotation matrix, we obtain a
% rotation that rotates vectors to the BPP based on the "reflected" B vector.
% If, for example, the original B vector resulted in a rotation angle (here, theta)
% of +120°, the adjusted rotation matrix, DMPA2BPP, describes the BPP coordinate
% system as rotated -60° from DSIz. All of this is just to minimize viewer
% discomfort interpreting the plot when B_DSIz changes sign.
%
% 2. Rotating one coordinate system in a direction (here, theta) relative to
% another coordinate system results in a requirement that elements in the original
% coordinatej system be rotated in the opposite direction to correctly specify
% their coordinates in the rotated system. The coordinates of a point or vector
% in the rotated coordinate system are now given by a rotation matrix which is
% the transpose of the coordinate axis rotation matrix.
%
% Additional notes:
% GDU locations must be rotated into BPP.
% The drift step and the target lie in BPP by definition, so they don't need to
% be rotated.
% Beams lie in planes parallel to BPP in BCS, and they all (ideally) point toward the
% drift step. When rotated into BPP, they remain in parallel planes.
% In this code segment, one plot shows BPP rotated relative to BCS, and another
% shows the GDU plane and beams rotated into BPP.

flipFlag = 1.0;
% if centerBeamB_u (3) < 0
% 	flipFlag = -1.0
% end

DMPA2BPPy = [ 0.0; 1.0; 0.0 ];
DMPA2BPPx = cross (DMPA2BPPy, centerBeamB_u);
DMPA2BPPx = DMPA2BPPx / norm (DMPA2BPPx, 2);
DMPA2BPPy = cross (centerBeamB_u, DMPA2BPPx);
DMPA2BPP  = [ DMPA2BPPx'; DMPA2BPPy'; centerBeamB_u' ];

% P = DMPA2BPP * p, where p is some point or vector in DMPA, can be interpreted as
% EITHER the fixed location of p expressed as P as seen in BPP
% OR p, expressed in DMPA coords, rotated in DMPA in the opposite magnitude and direction
% of the rotation of BPP from DMPA.
% This latter interpretation is the same as rotating P in BPP into DMPA coords.
% We see this here when we execute scBPPinDSI = DMPA2BPP * scBPP.

BPP2DMPA = DMPA2BPP';

% Plot the spacecraft outline in DMPA
if Use_OCS_plot
	set (0, 'CurrentFigure', hDMPA_plot) % hDMPA_plotElements
	hDMPA_plotElements (1) = plot3 ( ...
		instrPlaneDMPA (1,:), instrPlaneDMPA (2,:), instrPlaneDMPA (3,:), ...
		'LineStyle', '-', 'LineWidth', 1.0, 'Color', myDarkTeal);
	hold on
	% Plot the BPP plane in DMPA, along with the BPP normal vector (later).
	% Use the virtual sc outline as an arbitrary shape. REMEMBER, BPP is PERPENDICULAR
	% to B, and B is what we used to calculate the rotation matrix.
	BPPinDMPA = BPP2DMPA * instrPlaneDMPA;
	hDMPA_plotElements (2) = plot3 ( ...
		BPPinDMPA (1,:), BPPinDMPA (2,:), BPPinDMPA (3,:), ...
		'LineStyle', '-', 'LineWidth', 1.0, 'Color', 'blue');
end

set (0, 'CurrentFigure', hBPP_plot) % hBPP_plotElements DMPA2BPP
% Plot the sc as seen from BPP
instrPlaneDMPAinBPP = DMPA2BPP * instrPlaneDMPA;
% hBPP_plotElements (1) = plot3 (instrPlaneDMPAinBPP (1,:), instrPlaneDMPAinBPP (2,:), flipFlag * instrPlaneDMPAinBPP (3,:), 'LineStyle', '-', 'LineWidth', 1.0, 'Color', myDarkTeal);
hBPP_plotElements (1) = plot3 ( ...
	instrPlaneDMPAinBPP (1,:), ...
	instrPlaneDMPAinBPP (2,:), ...
	instrPlaneDMPAinBPP (3,:), ...
	'LineStyle', '-', 'LineWidth', 1.0, 'Color', myDarkTeal);
hold on
view ([   0 90 ]); % Azimuth, Elevation in degrees, std x-y plot

% Now plot the sc, which exists in DMPA, in BPP. The GDUs will appear on the
% periphery of this view of the sc from BPP. This may seem counter-intuitive,
% but imagine that the sc starts in BPP and rotates to its present position
% in DMPA... viewed in BPP vector space... that is why we use BBB2DSI here.
% This same line of reasoning is applied to the GDUs and the beams.
hBPP_plotElements (2) = plot3 (instrPlaneDMPA (1,:), instrPlaneDMPA (2,:), instrPlaneDMPA (3,:), 'LineStyle', '-', 'LineWidth', 1.0, 'Color', 'white');

if PlotHoldOff
	if (length (hBPP_plotElements) > 2)
		if Use_OCS_plot
			hDMPA_plotElements (3: length (hDMPA_plotElements)) = [];
		end
		hBPP_plotElements (3: length (hBPP_plotElements)) = [];
	end;
end

if Use_OCS_plot
	set (0, 'CurrentFigure', hDMPA_plot) % hDMPA_plotElements
	set (hDMPA_plot, 'color', 'white'); % sets the color to white
	AxisMax = 4;
	axis ([ -AxisMax AxisMax  -AxisMax AxisMax  -AxisMax AxisMax  -AxisMax AxisMax ]); % expanded axes for viewing larger drift steps
	axis square;
	axis vis3d;
	axis on;
	grid on;
	set (gca, 'XColor', myLightGrey4, 'YColor', myLightGrey4, 'ZColor', myLightGrey4)
end

set (0, 'CurrentFigure', hBPP_plot) % hBPP_plotElements DMPA2BPP
set (hBPP_plot, 'color', 'white'); % sets the color to white
AxisMax = 4;
axis ([ -AxisMax AxisMax  -AxisMax AxisMax  -AxisMax AxisMax  -AxisMax AxisMax ]); % expanded axes for viewing larger drift steps
axis square;
axis vis3d;
axis on;
grid on;
set (gca, 'XColor', myLightGrey4, 'YColor', myLightGrey4, 'ZColor', myLightGrey4)

EDP_dataInRange = true;
for BeamBracketIndex = (iiSorted_beam_t2k - beamsInWindow) : (iiSorted_beam_t2k + beamsInWindow)

	% BeamBracketIndex is uint32, so can't be < 0; just check for max size
	if ((BeamBracketIndex > 0) & (BeamBracketIndex < nBeams)) % catch array edge faults
	  iBeam = iSorted_beam_t2k (BeamBracketIndex);
		if (abs (centerBeam_t2k - edi_gd_beam_t2k (iBeam)) < (beamWindowSecs * 1e9)) % s ~> ns
			GDU_ID = edi_gd_ID (iBeam);
		  if iBeam == iCenterBeam
		    BeamColor = 'red';
		  else
		    BeamColor = myDarkGreen;
		  end

			GDU_Loc = edi_gd_virtual_dmpa (:, iBeam);
			if Use_OCS_plot
				set (0, 'CurrentFigure', hDMPA_plot) % hDMPA_plotElements
				if GDU_ID == 1
					hDMPA_plotElements (3) = plot3 (GDU_Loc (1), GDU_Loc (2), GDU_Loc (3), 'LineStyle', 'none', ...
						'Marker', 'o', 'MarkerFaceColor', 'white',   'MarkerEdgeColor', BeamColor, 'MarkerSize', 5.0); % GDU 1
				else
					hDMPA_plotElements (4) = plot3 (GDU_Loc (1), GDU_Loc (2), GDU_Loc (3), 'LineStyle', 'none', ...
						'Marker', 'o', 'MarkerFaceColor', myDarkGreen, 'MarkerEdgeColor', BeamColor, 'MarkerSize', 5.0); % GDU 2
				end
			end

			GDU_LocBPP = DMPA2BPP * GDU_Loc;
			set (0, 'CurrentFigure', hBPP_plot) % hBPP_plotElements DMPA2BPP
			disp (sprintf ('GDU: %d edi_gd_beam_t2k: %s x: %4.1f y: %4.1f', ...
				GDU_ID, datestr (edi_gd_beam_dn (iBeam), 'yyyy-mm-dd HH:MM:SS.fff'), ...
				GDU_LocBPP (1), GDU_LocBPP (2) ))
			if GDU_ID == 1
				hBPP_plotElements (3) = plot3 (GDU_LocBPP (1), GDU_LocBPP (2), GDU_LocBPP (3), 'LineStyle', 'none', ...
					'Marker', 'o', 'MarkerFaceColor', 'white',   'MarkerEdgeColor', BeamColor, 'MarkerSize', 5.0); % GDU 1
			else
				hBPP_plotElements (4) = plot3 (GDU_LocBPP (1), GDU_LocBPP (2), GDU_LocBPP (3), 'LineStyle', 'none', ...
					'Marker', 'o', 'MarkerFaceColor', myDarkGreen, 'MarkerEdgeColor', BeamColor, 'MarkerSize', 5.0); % GDU 2
			end

			% Technically, before we use this angle, we should adjust for drift speed angle change twixt outbound/inbound beam
			FiringDir = edi_gd_fv_dmpa (:, iBeam); % 3D

			BeamStart = GDU_Loc - 5.0 * FiringDir;
			BeamEnd   = GDU_Loc +       FiringDir;
			if Use_OCS_plot
				set (0, 'CurrentFigure', hDMPA_plot) % hDMPA_plotElements
				hDMPA_plotElements (5) = line ( [ BeamStart(1) BeamEnd(1) ], [ BeamStart(2) BeamEnd(2) ], [ BeamStart(3) BeamEnd(3) ], ...
					'LineStyle', ':', 'LineWidth', 1.0, 'Color', BeamColor); % Beam
			end

			BeamStartBPP =  DMPA2BPP * BeamStart;
			BeamEndBPP   =  DMPA2BPP * BeamEnd;
			set (0, 'CurrentFigure', hBPP_plot) % hBPP_plotElements BPP2DMPA
			% The beams are // already // parallel to BPP, so all we need to do is move them to the GDUs.
			% This should place them all in the same plane.
			hBPP_plotElements (5) = line ( [ BeamStartBPP(1) BeamEndBPP(1) ], [ BeamStartBPP(2) BeamEndBPP(2) ], [ BeamStartBPP(3) BeamEndBPP(3) ], ...
				'LineStyle', ':', 'LineWidth', 1.0, 'Color', BeamColor); % Beam

			% Find the most probable beam convergence for the drift step "target", S*
			% Theory:
			% Just above, the beams were rotated into BPP. Those beams are parallel to BPP (perpendicular to B),
			% though shifted in BBPz according to the tilt of the spacecraft in BPP.
			% If we assume a value of zero for the BPPz component of the beam, then we have but to solve
			% for the intersections of the beams in 2D. Here we gather the info that we'll need later:
			% slope and y-intercept. After processing all the beams, we'll calculate the intersections.
			% y = mx + b => m = dy/dx = (y2-y1) / (x2-x1); b = y1 - m*x1
% 			nTargetBeams = nTargetBeams + 1;
%% 			targetBeamsGDU_Loc (1:2, nTargetBeams) = GDU_LocBPP (1:2); % Effectively sets GDUz_loc to zero.
%% 			targetBeamAngle    (:, nTargetBeams)   = DMPA2BPP * FiringDir;
% 			targetBeam2DSlope (nTargetBeams) = (BeamEndBPP (2) - BeamStartBPP (2)) / (BeamEndBPP (1) - BeamStartBPP (1)); % dy/dx
% 			targetBeam2Db     (nTargetBeams) = BeamStartBPP (2) - targetBeam2DSlope (nTargetBeams) * BeamStartBPP (1); % b = y - mx
% 			FiringDirBPP (:, nTargetBeams)   = DMPA2BPP * FiringDir (:);

			% Find the EFW L2 target in the BPP plane defined by the centerBeamB and DMPA2BPP for this beam
			% 1. Rotate BavgSCS (nT) into BPP coordinates (needed for edp_E_interp calculations following)
			% 2. Rotate edp_E_interp (mV/m) into BPP coordinates
			% 3. Calculate the EFW BPP target location (meters) using (see inline references)
			% -------------------------------------------------------------------------------------

	 	  if iBeam == iCenterBeam
% disp 'iBeam == iCenterBeam'
				% centerBeamB : nT, 3D beam B-field vector, DMPA
				if Use_OCS_plot
					set (0, 'CurrentFigure', hDMPA_plot) % hDMPA_plotElements
					hDMPA_plotElements (6) = line ( [ 0.0 centerBeamB_u(1) ], [ 0.0 centerBeamB_u(2) ], [ 0.0 centerBeamB_u(3) ], ...
						'LineStyle', '-', 'LineWidth', 2.0, 'Color', myDarkBlue); % B field vector
				end

				set (0, 'CurrentFigure', hBPP_plot) % hBPP_plotElements DMPA2BPP
				% centerBeamB_u is the real B unit vector.
				CenterBeamBavgBPP_u = DMPA2BPP * centerBeamB_u;
				hBPP_plotElements (6) = line ( [ 0.0 CenterBeamBavgBPP_u(1) ], [ 0.0 CenterBeamBavgBPP_u(2) ], [ 0.0 CenterBeamBavgBPP_u(3) ], ...
					'LineStyle', '-', 'LineWidth', 2.0, 'Color', myDarkBlue); % B field vector

				if Use_OCS_plot
					set (0, 'CurrentFigure', hDMPA_plot) % hDMPA_plotElements
					hDMPA_plotElements (7) = line ( ...
						[ 0.0 centerBeamEu(1) ], [ 0.0 centerBeamEu(2) ], [ 0.0 centerBeamEu(3) ], ...
						'LineStyle', '-', 'LineWidth', 2.0, 'Color', myDarkRed); % E DMPA field unit vector
				end

				set (0, 'CurrentFigure', hBPP_plot) % hBPP_plotElements DMPA2BPP
				centerBeamEu_bpp = DMPA2BPP * centerBeamEu;
				hBPP_plotElements (7) = line ( ...
					[ 0.0 centerBeamEu_bpp(1) ], [ 0.0 centerBeamEu_bpp(2) ], [ 0.0 centerBeamEu_bpp(3) ], ...
					'LineStyle', '-', 'LineWidth', 2.0, 'Color', myDarkRed); % E DMPA field vector

				if S_star_dmpa (1) ~= -edi_d_dmpa_fillVal % S* is the negative of the drift step, so FILLVALs are negated, too
					if Use_OCS_plot
						set (0, 'CurrentFigure', hDMPA_plot) % hDMPA_plotElements
						hDMPA_plotElements (8) = plot3 ( ...
							S_star_dmpa (1), S_star_dmpa (2), S_star_dmpa (3), ...
							'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red', 'MarkerSize', 5.0);
					end

					set (0, 'CurrentFigure', hBPP_plot) % hBPP_plotElements
					S_star_bpp = DMPA2BPP * S_star_dmpa;
					hBPP_plotElements (8) = plot3 ( ...
						S_star_bpp (1), S_star_bpp (2), S_star_bpp (3), ...
						'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'red', 'MarkerEdgeColor', 'red', 'MarkerSize', 5.0);
				end

				% Perform E·B = 0 in 3D to get the full EFW perp vector in DMPA (need z-axis)
				% Find 2 EFW records whose ssm brackets the center beam; interpolate the E-fields wrt ssm
				iEDP_t2kHI = find ((edp_t2k - edi_gd_beam_t2k (iBeam)) > 0.0, 1); % First edp_t2k > center beam time

				% EFW L2 ssm must be close enough to beam time to be useful
				% We don't want to interpolate over more than 0.1 s.
				% 0.1s / EDP_samplePeriod (0.0078125s @ 128 Hz) ~> 12.8 records
				% tt2000 is in ns: 0.1s ~> 1e8 ns. Because there are potentially
				% ~13 records in this amount of time, there is no loss of generality
				% seting the limit at 0.1s.
				if abs ((edp_t2k (iEDP_t2kHI) - edi_gd_beam_t2k (iBeam)) < 1.0e8)
					iEDP_t2kLO = iEDP_t2kHI - 1; % First edp_t2k < center beam time
					if iEDP_t2kLO == 0 % We must bump this if we are looking at the first record... maybe we should start @ 2?
						iEDP_t2kLO = 1;
						iEDP_t2kHI = 2;
					end

					edp_t2kLO = edp_t2k (iEDP_t2kLO);
					disp (sprintf ('EDP LO time: %s', char(spdfencodett2000(edp_t2kLO)) )) % V&V

					edp_t2kHI = edp_t2k (iEDP_t2kHI);
					disp (sprintf ('EDP HI time: %s', char(spdfencodett2000(edp_t2kHI)) )) % V&V

					% Remember that there are ~13 recs / 0.1s @ 128 Hz sample rate
					if (edp_t2kHI - edp_t2kLO) > 1.0e8
						warning ('Interpolating EFW data points more than 100 ms apart.');
					end

					interp_frac  = double (centerBeam_t2k - edp_t2kLO) / double (edp_t2kHI - edp_t2kLO);

% Don't remove until decision re E.B=0
% 					edp_ExLO = edp_E3D_dsl (1, iEDP_t2kLO);
% 					edp_EyLO = edp_E3D_dsl (2, iEDP_t2kLO);
% 					edp_EzLO = - (edp_ExLO * centerBeamB (1) + edp_EyLO * centerBeamB (2)) / centerBeamB (3);
% 					edp_ELO  = [ edp_ExLO; edp_EyLO; edp_EzLO ];

					edp_ELO  = edp_E3D_dsl (:, iEDP_t2kLO);
	% 				disp ( [ '....... E-field LO: ', sprintf('%.3f %.3f %.3f',edp_ELO) ] ); % V&V

% 					edp_ExHI = edp_E3D_dsl (2, iEDP_t2kHI);
% 					edp_EyHI = -edp_E3D_dsl (1, iEDP_t2kHI);
% 					edp_EzHI = - (edp_ExHI * centerBeamB (1) + edp_EyHI * centerBeamB (2)) / centerBeamB (3);
% 					edp_EHI  = [ edp_ExHI; edp_EyHI; edp_EzHI ];
					edp_EHI = edp_E3D_dsl (:, iEDP_t2kHI);
	% 				disp ( [ '....... E-field HI: ', sprintf('%.3f %.3f %.3f',edp_EHI) ] ); % V&V

					edp_E_diff = edp_EHI - edp_ELO;
	% 				disp ( [ '....... E-field HI-LO delta: ', sprintf('%.3f %.3f %.3f',edp_E_diff) ] ); % V&V

					edp_E_interp   = edp_ELO + edp_E_diff * interp_frac; % The interpolated value of edp_E at EDI beam time
					edp_E_interp_u = edp_E_interp / norm (edp_E_interp, 2);
					disp ( [ 'edp_E_interp DSL : ', sprintf('%+8.3f %+8.3f %+8.3f', edp_E_interp) ] );
					disp ( [ 'centerBeamB  DMPA: ', sprintf('%+8.3f %+8.3f %+8.3f', centerBeamB) ] );

					if Use_OCS_plot
						set (0, 'CurrentFigure', hDMPA_plot) % hDMPA_plotElements
						hDMPA_plotElements (9) = line ( [ 0.0 edp_E_interp_u(1) ], [ 0.0 edp_E_interp_u(2) ], [ 0.0 edp_E_interp_u(3) ], ...
							'LineStyle', '-', 'LineWidth', 2.0, 'Color', myDarkGreen); % E DMPA field vector
					end

					set (0, 'CurrentFigure', hBPP_plot) % hBPP_plotElements DMPA2BPP
					edp_E_BPPu = DMPA2BPP * edp_E_interp_u;
					hBPP_plotElements (9) = line ( [ 0.0 edp_E_BPPu(1) ], [ 0.0 edp_E_BPPu(2) ], [ 0.0 edp_E_BPPu(3) ], ...
						'LineStyle', '-', 'LineWidth', 2.0, 'Color', myDarkGreen); % E DMPA field vector

					% [9] Chen, [2-3]
					% [12] Messeder, (validation calculations), 2011
					gyroFrequency = (-q * centerBeamB2n * nT2T) / e_mass; % (SI) (|q| is positive here.)
					gyroPeriod = (2.0 * pi / gyroFrequency);      % (SI) The result is usually on the order of a few ms

					% We have filtered for beam loops = 1 (N =1)
					% d = (E x B / B^2) * N * beam period ???
					% [9] Chen, [2-15]
					% [1] Paschmann, hDMPA_plotElements 241, (3)
					% [1] p244: The drift vector is directed from the target (S*) to the detector.
					% ExB_targetBPP = SCS2BPP * ( cross (edp_E_interp * mV2V, BavgSCS * nT2T) / (centerBeamB2n * nT2T)^2 ) * gyroPeriod % (meters) with SI conversions

					% vE = v in direction of E; d = drift step; T = gyroPeriod
					% ( vE = d/T ) = ExB/|B|^2 ~> d / T * |B|^2 = ExB --- Pacshmann, 1998, 2001, EDI for Cluster
					% Cross from the left with B: B x [] = BxExB
					% where BxExB = E(B dot B) - B(E dot B)
					% The second term is zero because we are assuming E is perpendicular to B.
					% B x [ d/T * |B|^2 = E * |B|^2 ~> E = B x d/T

					% The virtual source S* is the NEGATIVE of the real drift step; S* is an imaginary point
					% This should be E = B x v, but B, v are swapped here because we need the real drift step (drift velocity),
					% not the virtual source, S*. See relevant publications on Cluster drift step
					% and 'EDI_beams_and_virtual_source_demo_0101.m'.

% 					edp_E_interp_d_dmpa = ((cross (edp_E_interp, centerBeamB) / centerBeamB2n^2 ) * C2V2T * gyroPeriod); % (m), simplified
% 					edp_E_interp_d_dmpa_u = edp_E_interp_d_dmpa / norm (edp_E_interp_d_dmpa, 2);

% !!! what should this be?
% 					edp_E_interp_d_bpp = DMPA2BPP * edp_E_interp_d_dmpa;
% 					if (abs (edp_E_interp_d_bpp (3)) > 1.0e-4)  % we really need only this one...
% 						warning ('|edp_E_interp_d_bpp_z| > 1.0e-4')
% 					end

					edp_zoomStart = iEDP_t2kLO - 100;
					edp_zoomEnd   = edp_zoomStart + 200;
					MEEdrift_plot_EDP_zoomed_region

				else
					disp 'Warning!! Interpolating EFW data points more than 100 ms apart.'
					EDP_dataInRange = false;
				end

% 				if S_star_dmpa (1) ~= -edi_d_dmpa_fillVal % S* is the negative of the drift step, so FILLVALs are negated, too
% 					plot3 (S_star_dmpa (1), S_star_dmpa (2), S_star_dmpa (3), ...
% 						'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'green', 'MarkerSize', 10.0);
% 					disp ( [ 'S* BPP          : ', sprintf('%+8.3f %+8.3f %+8.3f', S_star_bpp) ] );
% 				end
			end % if BeamBracketIndex == iCenterBeam
		end % if beams within beamWindowSecs seconds of center beam

	end % if ((BeamBracketIndex > 0) & (BeamBracketIndex < BeamRecords))
end % BeamBracketIndex

%{
% Find center of beam intersections, if there are > 1 beams
% There is a chance that some lines will be parallel, or nearly so... this must be eventually handled
% if this method proves viable. y = mx + b => -mx + y = b
% | -m 1| |x| = |b|, where || signals a matrix; then [x y] = M \ b
set (0, 'CurrentFigure', hBPP_plot) % hBPP_plotElements
if nTargetBeams > 1
	nbeamIntercepts = 0;
	for i = 2: nTargetBeams
		for j = 1: i-1
			XY = [ ...
				-targetBeam2DSlope(i) 1 ;
				-targetBeam2DSlope(j) 1 ];
			b = [ targetBeam2Db(i) ; targetBeam2Db(j) ];
			nbeamIntercepts = nbeamIntercepts + 1;
			beamIntercepts (:, nbeamIntercepts) = XY \ b;
		end
	end
	beamIntercepts (:, nbeamIntercepts+1:end) = [];
	beamInterceptMean = mean (beamIntercepts, 2)
	beamInterceptStdDev = std (beamIntercepts, 1, 2)

% 	ibx = find ( abs(beamIntercepts (1, :) - beamInterceptMean (1)) > 2.0 * beamInterceptMean (1) );
% 	beamIntercepts (:,ibx) = [];
% 	ibx = find ( abs(beamIntercepts (2, :) - beamInterceptMean (2)) > 2.0 * beamInterceptMean (2) );
% 	beamIntercepts (:, ibx) = [];

	nBeamIntercepts = length (beamIntercepts);
	if nBeamIntercepts > 3
		ibx = find ( abs(beamIntercepts (1, :) - beamInterceptMean (1)) > 1.0 * beamInterceptStdDev (1) ); % n-sigma acceptance
		beamIntercepts (:,ibx) = [];
		ibx = find ( abs(beamIntercepts (2, :) - beamInterceptMean (2)) > 1.0 * beamInterceptStdDev (2) );
		beamIntercepts (:, ibx) = [];

		nBeamIntercepts = length (beamIntercepts);
		if nBeamIntercepts > 3
			beamInterceptMean = mean (beamIntercepts, 2)
			beamInterceptStdDev = std (beamIntercepts, 1, 2)
% 			ibx = find ( abs(beamIntercepts (1, :) - beamInterceptMean (1)) > 2.0 * beamInterceptStdDev (1) );
% 			beamIntercepts (:,ibx) = [];
% 			ibx = find ( abs(beamIntercepts (2, :) - beamInterceptMean (2)) > 2.0 * beamInterceptStdDev (2) );
% 			beamIntercepts (:, ibx) = [];
		end
% 		beamInterceptMean = mean (beamIntercepts, 2)
% 		beamInterceptStdDev = std (beamIntercepts, 1, 2)
	end

	if plotBeamDots
		plot3 (beamIntercepts (1,:), beamIntercepts (2,:), 0.0*[1:nBeamIntercepts], ...
			'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', myLightGrey4, 'MarkerEdgeColor', myLightGrey4, 'MarkerSize', 2.0);
	end

% 	plot3 (beamInterceptMean (1), beamInterceptMean (2), 0.0, ...
% 		'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'green', 'MarkerSize', 10.0);
% 	bweighted = mean (beamIntercepts, 2)

	plot3 (beamInterceptMean (1), beamInterceptMean (2), 0.0, ...
		'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'blue', 'MarkerEdgeColor', 'green', 'MarkerSize', 10.0);
	disp ( [ 'Beam convergence: ', sprintf('%+8.3f %+8.3f %+8.3f', beamInterceptMean, 0.0) ] );
end
%}

if Use_OCS_plot
	if (length (hDMPA_plotElements) > 9) hDMPA_plotElements (9:length (hDMPA_plotElements)) = []; end;
end
if (length (hBPP_plotElements) > 9) hBPP_plotElements (9:length (hBPP_plotElements)) = []; end;

if Use_OCS_plot
	set (0, 'CurrentFigure', hDMPA_plot) % hDMPA_plotElements
	hDMPA__BPP_legendAxes = gca;
	p = hDMPA_plotElements;
	PlotView = '3D'; % flag for ...annotate...
	% 	MEEdrift_annotate_DMPA_BPP_plots
	% xlabel ('x')
	% ylabel ('y')
	figure (hDMPA_plot);
	view ([ 115 20 ]); % Azimuth, Elevation in degrees, 3D view
	zoom (0.8)
end

set (0, 'CurrentFigure', hBPP_plot) % hBPP_plotElements DMPA2BPP
hDMPA__BPP_legendAxes = gca;
p = hBPP_plotElements;
PlotView = '2D';
% 	MEEdrift_annotate_DMPA_BPP_plots
% xlabel ('x')
% ylabel ('y')

figure (hBPP_plot);
view ([   0 90 ]); % Azimuth, Elevation in degrees, std x-y plot

if EDP_dataInRange
	if Use_OCS_plot
		set (0, 'CurrentFigure', hDMPA_plot) % hDMPA_plotElements
		hDMPA__BPP_legendAxes = gca;
		p = hDMPA_plotElements;
		PlotView = '3D'; % flag for ...annotate...
		MEEdrift_annotate_DMPA_BPP_plots

		figure (hDMPA_plot);
		view ([ 115 20 ]); % Azimuth, Elevation in degrees, 3D view
		zoom (0.8)
	end

	set (0, 'CurrentFigure', hBPP_plot) % hBPP_plotElements DMPA2BPP
	hDMPA__BPP_legendAxes = gca;
	p = hBPP_plotElements;
	PlotView = '2D';
	MEEdrift_annotate_DMPA_BPP_plots

	figure (hBPP_plot);
	view ([   0 90 ]); % Azimuth, Elevation in degrees, std x-y plot
% 	TightFig;

	% Update the EFW data plot index line
	figure (hEDP_plot);
	axes (hEDP_plot_mainAxes);
	delete (hEDP_plot_index_line);
	EDP_dataPlotIndexLine = centerBeam_dn;
disp ( sprintf ('Skip via EDP: EDP_dataPlotIndexLine %s', datestr (EDP_dataPlotIndexLine, 'yyyy-mm-dd HH:MM:ss.fff') ) ) % V&V
% 	datestr (centerBeam_dn, 'yyyy-mm-dd HH:MM:ss.fff') % V&V
	hold on
	hEDP_plot_index_line = line ( [EDP_dataPlotIndexLine EDP_dataPlotIndexLine], get (hEDP_plot_mainAxes, 'YLim'), 'Color', 'red' , 'LineStyle', '-' , 'LineWidth', 2);

	% update EFWplotMagnifier here
% 	iLeftMag_ssm  = max (1,                  find (edp_t2k >= (BdvE_dn - 4.0), 1) );
% 	iRightMag_ssm = min (nEDP, find (edp_t2k >= (BdvE_dn + 4.0), 1) );

% 	MMS_TimeMagSeries {1}  = (edp_t2k (1, iLeftMag_ssm: iRightMag_ssm)) / 86400.0; % make ssm (s) > DateNum (fractional days)
% 	MMS_DataMagSeries {1}  = edp_E3D_dsl  (:, iLeftMag_ssm: iRightMag_ssm);

% 	if FilterData
% 		MMS_DataMagSeriesf {1} = EFW_L2_E2Df (:, iLeftMag_ssm: iRightMag_ssm);
% 	else
% 		MMS_DataMagSeriesf {1} = edp_E3D_dsl (:, iLeftMag_ssm: iRightMag_ssm);
% 	end

% 	axes (hEFWsubplotMagAxes);
% 	delete (hEFWplotMag);
% 	hEFWplotMag = plot (hEFWsubplotMagAxes, MMS_TimeMagSeries {1}, MMS_DataMagSeries {1}, MMS_TimeMagSeries {1}, MMS_DataMagSeriesf {1}, 'LineStyle', '-', 'Marker', '.', 'MarkerSize', 2);

% 	DateStartMag = MMS_TimeMagSeries {1} (1);
% 	DateEndMag   = MMS_TimeMagSeries {1} (end);
% 	xlim ([DateStartMag DateEndMag]);
% 	xTickData = linspace (DateStartMag, DateEndMag, 5);
% 	set (gca, 'XTick', xTickData);
% 	datetick ('x', 'HH:MM:SS','keepticks', 'keeplimits');
% 	ZoomDateTicks ('on');
%}
else
  disp 'Warning!!! No EDP coordinated data in range!!!'
end
