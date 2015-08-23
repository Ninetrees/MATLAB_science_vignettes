function [] = edi_drift_step_plot ( ...
	obsID, ...
	B_t2k, ...
	gd_ID, ...
	gd_virtual_bpp, ...
	gd_fv_bpp, ...
	B_sdcs, ...
	S_star_bpp, ...
	SDCS2BPP, ...
	GrubbsBeamIntercepts, GrubbsBeamInterceptMean, GrubbsBeamInterceptMean_stdDev, ...
	P0, dsQuality);

	myLibAppConstants % custom colors; set default axis colors
	cPathSep = pathsep;
	cFileSep = filesep;
	global dotVersion pause_each_plot plot_dots save_beam_plots show_unit_vectors
	persistent hBPP_figure; isempty (hBPP_figure); % seems to be needed to trigger proper logic below
% 	disp 'entering edi_drift_step_plot'

	if isempty (hBPP_figure)
		hBPP_figure = figure ('Position', [ 400   150   850   800 ]); set (hBPP_figure, 'Visible', 'off')
		set (hBPP_figure, 'WindowStyle', 'normal')
		set (hBPP_figure, 'DockControls', 'off')
		set (gcf, 'name', 'B_E_field_GDU_locs_beams', 'NumberTitle', 'off', 'visible', 'on');
		set (0, 'CurrentFigure', hBPP_figure) % hBPP_plotElements
	end

	% see also clf (fig, 'reset');
	if pause_each_plot % next commands work best when on the same line !!!
		clf (hBPP_figure); set (hBPP_figure, 'Visible', 'on');
	else
		clf (hBPP_figure); set (hBPP_figure, 'Visible', 'off');
	end

	EDI1gunLoc = [ -1.45598,  1.11837, 0.0 ]; % EDI2 gun atan2(-1.11837, 1.45598)*180/pi ~> -37.52865°
	EDI1detLoc = [ -1.35885,  1.03395, 0.0 ]; % EDI2 det atan2(-1.03395, 1.35885)*180/pi ~> -37.26753°
	EDI2gunLoc = [  1.45598, -1.11837, 0.0 ]; % EDI1detLoc:EDI1gunLoc angle = atan2(1.11837-1.03395, -1.45598+1.35885)*180/pi ~> -40.995°
	EDI2detLoc = [  1.35885, -1.03395, 0.0 ]; % norm(EDI1gunLoc-EDI1detLoc,2) = 0.128689
	mmsEDI_VirtualRadius = norm (EDI1gunLoc-EDI2detLoc, 2);

	theta = 0.0: 1.0: 360.0;
	virtual_instr_plane = mmsEDI_VirtualRadius * ...
		[ cosd(theta); sind(theta); zeros(1,length(theta),'double') ];

	BPP_plane = virtual_instr_plane; % re-use v_i_p, because the BPP_plane looks normal in BPP
	hBPP_plotElements (1) = plot3 ( ...
		BPP_plane (1,:), ...
		BPP_plane (2,:), ...
		BPP_plane (3,:), ...
		'LineStyle', '-', 'LineWidth', 1.0, 'Color', myDarkBlue);

	AxisMax = 4;
	axis ([ -AxisMax AxisMax  -AxisMax AxisMax  -AxisMax AxisMax ]); % expanded axes for viewing larger drift steps
	axis square
	axis vis3d
	axis on
	grid on
	set (gca, 'XColor', myLightGrey4, 'YColor', myLightGrey4, 'ZColor', myLightGrey4)
	hold on % must not be set until AFTER first item is plotted.
	xlabel ('X (m)');
	ylabel ('Y (m)');
	zlabel ('Z');
	legend ('-DynamicLegend'); % no effect on legend problem
	view ([ 0 90 ])
	set (gcf, 'Units', 'normal')
	set (gca, 'Position', [0.01 0.01 0.99 0.99])

	GDU_planeInBPP = SDCS2BPP * virtual_instr_plane; % instrument plane rotated into BPP
	hBPP_plotElements (2) = plot3 ( ...
		GDU_planeInBPP (1,:), ...
		GDU_planeInBPP (2,:), ...
		GDU_planeInBPP (3,:), ...
		'LineStyle', '-', 'LineWidth', 1.0, 'Color', myDarkRed);
	% disp 'step 2 ~> DMPA in red'

	igd_ID = find (gd_ID == 1); % GDU 1
	if ~isempty (igd_ID)
		p = plot3 ( ...
			gd_virtual_bpp (1,igd_ID), ...
			gd_virtual_bpp (2,igd_ID), ...
			gd_virtual_bpp (3,igd_ID), ...
			'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'white',   'MarkerEdgeColor', myDarkGreen, 'MarkerSize', 5.0);
	else % plot nothing for legend
		p = plot3 ( 0.0, 0.0, 0.0, ...
			'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'white',   'MarkerEdgeColor', 'white', 'MarkerSize', 5.0);
	end
	hBPP_plotElements (3) = p;

	igd_ID = find (gd_ID == 2); % GDU 2
	if ~isempty (igd_ID)
		p = plot3 ( ...
			gd_virtual_bpp (1,igd_ID), ...
			gd_virtual_bpp (2,igd_ID), ...
			gd_virtual_bpp (3,igd_ID), ...
			'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', myDarkGreen, 'MarkerEdgeColor', myDarkGreen, 'MarkerSize', 5.0);
	else % plot nothing for legend
		p = plot3 ( 0.0, 0.0, 0.0, ...
			'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'white',   'MarkerEdgeColor', 'white', 'MarkerSize', 5.0);
	end
	hBPP_plotElements (4) = p;
	% disp 'step 3 ~> GDUs as seen in BPP'

	for i=1: size (gd_virtual_bpp, 2)
		% The beams are // already [supposed to be] // parallel to BPP, so all we need to do is move them to the GDUs.
		% This should place them all in the same plane.
		GDU_Loc      = gd_virtual_bpp (:,i);
		BeamStartBPP = GDU_Loc - 6.0 * gd_fv_bpp (:,i);
		BeamEndBPP   = GDU_Loc + 6.0 * gd_fv_bpp (:,i);
% 		BeamStartBPP = BeamStart; % OCS2BPP * BeamStart;
% 		BeamEndBPP   = BeamEnd; % OCS2BPP * BeamEnd;

		p = plot3 ( ...
			[ BeamStartBPP(1) BeamEndBPP(1) ], ...
			[ BeamStartBPP(2) BeamEndBPP(2) ], ...
			[ BeamStartBPP(3) BeamEndBPP(3) ], ...
			'LineStyle', ':', 'LineWidth', 1.0, 'Color', myDarkBlue); % Beam
	end
	hBPP_plotElements (5) = p; % this doesn't fix the legend problem
	% disp 'step 4 ~> firing vectors as seen in BPP'

	% The virtual source S* and confidence box
	zeroVector = [ 0.0 ];
	hBPP_plotElements (6) = plot3 (GrubbsBeamInterceptMean (1), GrubbsBeamInterceptMean (2), zeroVector (1), ...
		'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', myOrange, 'MarkerEdgeColor', myOrange, 'MarkerSize', 5.0);

	set (hBPP_plotElements (6), 'Color', myOrange);

% 	disp ( sprintf ('Beam convergence: %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f %+8.3f', ...
% 		GrubbsBeamInterceptMean, GrubbsBeamInterceptStdDev, GrubbsBeamInterceptMean_stdDev) )

	ZConfidenceBounds = norminv (1.0 - (1.0 - P0) / 2.0);
	ConfidenceIntervals = [ ...
		GrubbsBeamInterceptMean-(ZConfidenceBounds*GrubbsBeamInterceptMean_stdDev),... % x, y lower limits
		GrubbsBeamInterceptMean+(ZConfidenceBounds*GrubbsBeamInterceptMean_stdDev) ];  % x, y upper limits
	disp (['Grubbs 84% confidence intervals   = ', sprintf('( %g, %g )', ConfidenceIntervals) ])

	if ( ((ConfidenceIntervals(1,2) - ConfidenceIntervals(1,1)) > 0.0) & ...
	     ((ConfidenceIntervals(2,2) - ConfidenceIntervals(2,1)) > 0.0) )
		r = rectangle ('Position', [ ...
			ConfidenceIntervals(1,1), ...
			ConfidenceIntervals(2,1), ...
			ConfidenceIntervals(1,2)-ConfidenceIntervals(1,1), ...
			ConfidenceIntervals(2,2)-ConfidenceIntervals(2,1) ], ...
			'LineStyle', '--', 'LineWidth', 1);
		set (r, 'edgecolor', myOrange)
	end

	S_star_bpp = [ GrubbsBeamInterceptMean; 0.0 ];
	if show_unit_vectors
		% B unit vector
		B_bpp = SDCS2BPP * B_sdcs;
		B_bpp_u = B_bpp / norm (B_bpp, 2);
		hBPP_plotElements (7) = line ( [ 0.0 B_bpp_u(1) ], [ 0.0 B_bpp_u(2) ], [ 0.0 B_bpp_u(3) ], ...
			'LineStyle', '-', 'LineWidth', 2.0, 'Color', MMS_plotColorx); % B field vector

		virtualSource_bpp_u = S_star_bpp / norm (S_star_bpp, 2);
		hBPP_plotElements (8) = line ( ...
			[ 0.0 virtualSource_bpp_u(1) ], ...
			[ 0.0 virtualSource_bpp_u(2) ], ...
			[ 0.0 virtualSource_bpp_u(3) ], ...
			'LineStyle', '-', 'LineWidth', 2.0, 'Color', MMS_plotColory); % B field vector
	end

	if plot_dots
		% Note here that NOT ALL intersections are plotted,
		% just those returned by Grubbs...()
% 		GrubbsBeamInterceptsSize = size (GrubbsBeamIntercepts)
% 		size (GrubbsBeamIntercepts, 2)
		p = scatter ( ...
			GrubbsBeamIntercepts (1,:), ...
			GrubbsBeamIntercepts (2,:), ...
			9.0, 'r', '*');
% 			5.0, myGold, '*');
% 			'MarkerFaceColor', myGold, 'MarkerEdgeColor', myGold);
% 		for i = 1: size (GrubbsBeamIntercepts, 2)
% 	 		fprintf ('.');
% 			p = plot3 ( ...
% 				GrubbsBeamIntercepts (1,1), ...
% 				GrubbsBeamIntercepts (2,1), ...
% 				0.0, ...
% 				'LineStyle', 'none', 'Marker', '*', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 2.0);
% 				'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', myGold, 'MarkerEdgeColor', myGold, 'MarkerSize', 2.0);
% 		end
% size (p)

		hBPP_plotElements (9) = p;
	end

	nBeamIntercepts = size (GrubbsBeamIntercepts, 2);
	bTimeStr = datestr (spdftt2000todatenum (B_t2k), 'yyyy-mm-dd HH:MM:ss');
	title ( { ...
		[ 'EDI drift step using beam convergence, ', dotVersion, ', B_{AVG} BPP'];
		[ 'MMS', obsID, ': ', bTimeStr, ', P_{0} = ', num2str(P0, 3), ', Points = ', num2str(nBeamIntercepts), ', Q = ', num2str(dsQuality) ];
		[ 'B_{DMPA} [',...
		  num2str(round (B_sdcs (1))),',',...
		  num2str(round (B_sdcs (2))),',',...
		  num2str(round (B_sdcs (3))), '] nT,   ', ...
		  'd_{BPP} [',...
		  num2str(S_star_bpp (1),'%5.2f'),',',...
		  num2str(S_star_bpp (2),'%5.2f'),',',...
		  num2str(S_star_bpp (3),'%5.2f'), '] m', ...
		]; ...
		[ '~> click plot to advance...'];
		[] }, 'Fontname', 'Times');

	virtualSource_bpp_angle = atan2d (S_star_bpp(2), S_star_bpp(1));
	if ( (virtualSource_bpp_angle > 0) & (virtualSource_bpp_angle < 90) ) % target in NorthEast
		legend_location = 'NorthWest';
	else
		legend_location = 'NorthEast';
	end
	switch length (hBPP_plotElements) % !!! Specifying legend (hBPP_plotElements, ... FINALLY fixed problem w legend not displaying correctly
		case 6
			legend (hBPP_plotElements, 'BPP plane', 'GDU plane in BPP', 'GDU 1', 'GDU 2', 'EDI beam', 'Virtual Source', ...
				'Location', legend_location );
		case 8
			legend (hBPP_plotElements, 'BPP plane', 'GDU plane in BPP', 'GDU 1', 'GDU 2', 'EDI beam', 'Virtual Source', ...
				'B_{BPP} (unit)', 'S*_{BPP} (unit)', ...
				'Location', legend_location );
		case 9
			legend (hBPP_plotElements, 'BPP plane', 'GDU plane in BPP', 'GDU 1', 'GDU 2', 'EDI beam', 'Virtual Source', ...
				'B_{BPP} (unit)', 'S*_{BPP} (unit)', 'Intercepts', ...
				'Location', legend_location );
	end
	if save_beam_plots
		% Example: 'mms2_edi_20150509_drift_step_BPP_163130_v1.03a.png'
		SavePlotFilename = [ ...
			'.' cFileSep 'mms', obsID, '_edi_', ...
			bTimeStr(1:4) bTimeStr(6:7) bTimeStr(9:10) '_drift_step_BPP', ...
			bTimeStr(12:13) bTimeStr(15:16) bTimeStr(18:19), ...
			'_', dotVersion, ...
			'_i', num2str(size(GrubbsBeamIntercepts, 2), '%04d'), '.png', ...
		];
		hgexport (gcf, SavePlotFilename, EDI_presentation_beam_plot_style);
	end

	if pause_each_plot
		dummy = waitforbuttonpress;
	end
end

% for i=1:3
% 	disp (sprintf (' i = %d\n', i))
% 	if i > 2
% 		disp 'i>2'
% 	else
% 		disp 'i <= 2'
% 		if i > 1
% 			disp 'i>1'
% 		else
% 			disp 'i=1'
% 		end
% 	end
% end
% disp 'end'
