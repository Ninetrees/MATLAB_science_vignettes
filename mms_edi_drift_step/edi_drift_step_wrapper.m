% edi_drift_step
%
% Purpose
% Using data from MMS EDI and mag gear, determine the drift step in BPP, the
% uncertainty in the drift step, the drift velocity, and the E-field derived from
% the drift step, using beam convergence statistics.
%
% This program is a test bed for production software. Routines vetted here are
% embedded in MMS production software.
%
% Calling Sequence(s):
% edi_drift_step_wrapper ()
% There are no parms passed in. Each version of the program reads specific MMS
% data files that may change from version to version.
% Data files for this version:
% 'mms2_edi_slow_ql_efield_20150509_v0.1.4.cdf'
% 'mms2_edp_comm_ql_dce2d_20150509000000_v0.1.0.cdf'
%
% Output:
% Optionally produces a series of drift step plots.
% Optionally saves drift step plots.
% Optionally produces a summary plot of B, SDP, and drift step E-field.
% ~~~> global program control flag notes
% pause_each_plot   ~ During plotting of beam sets, pauses until user clicks on plot.
%                   ~ Requires plot_beams
%                   ~ Plots are not visible if they are not paused
% plot_beams        ~ Causes plotting of 5s beam sets. If this is not set, drift steps
%                     are calculated, but not plotted.
% plot_dots         ~ Plot dots at beam intersections that are used in final drift step
%                     calculations.
%                   ~ Requires plot_beams
% plot_edi_E_sdcs   ~ Plot E-field from CDF file to compare with edi_drift_step (inactive).
% plot_summary      ~ Plot the summary plot (B, EDP E, drift step E, Quality, Quality histo).
% plot_sinx_weights ~ Plot the weighting function and quality crossover points.
% save_beam_plots   ~ Saves each 5s beam set to files.
%                     Ex: mms2_edi_20150509_drift_step_BPP163140_v1.05.00_i0001.png
%                   ~ Requires plot_beams
% use_v10502        ~ Remove parallel beams; otherwise, leave them in, but assign
%                     low weights to them.
%
% Important !!! notes about variable naming conventions
% iVarName is an index to VarName
% iiVarName is an index to an index to VarName; e.g., VarName (iVarName (iiVarName))
% nVarName is the number of elements in VarName
% _dn indicates a MATLAB datenum variable
% _frunc indicates fractional uncertainty: SD / mean
% _sdcs <~~~~<< is 'some despun coord system'. It stands for any that is useful.
% _t2k indicates a CDF TT2000 variable
% _u or ...u is a unit vector; _u used for clarity when needed
% 2n indicates a 2-norm result.
% ISA indicates intersection angle
% CI - confidence interval
% SD - std dev
% SDOM - std dev of the mean. See important notes in 'Accuracy, Error, Precision,
%  and Uncertainty.txt' regarding the use of the standard deviation (here, SD)
%  and the standard deviation of the mean (here, SDOM).
%  See http://onlinestatbook.com/2/estimation/mean.html ~> mu +- Z.nn * SDOM
% Percentile to Z-Score for a normal distribution
%  z score, also z critical value. See 'Accuracy, etc.'
%
% MATLAB release(s) MATLAB 8.3.0.532 (R2014a)
% Required Products None
%
% Plots that prompt questions:
%
% History:
% 2015-10- ~ v02.00.00:
% ~ Add uncertainty notes; lay the foundation for uncertainty propagation.
% ~ Standardize on SD and SDOM.
% ~ Clean up interface twixt edi_drift_step and edi_drift_step_plot.
% ~ edi_drift_step: develop better flow of var names for S_star and related vars.
% ~ Research statistics terms to include here; include references.
% ~ Update 'Accuracy, Error, Precision, and Uncertainty.txt'.
% ~ Propagate S* and B errors through to E-field results.
%
% 2015-09-04 ~ v01.06.00:
% ~ Important !!! Notes about variable naming conventions, UPDATED.
% ~ Change local, not CDF, vars _dmpa to _sdcs <~~~~<<.
% ~ Set minimum number of beams to use at 2..
% ~ Correct EDP data call to use obsID.
% ~ Change large plot types from lines to dots.
%
% 2015-08-11 ~ v01.05.06:
% ~ minor internal docs
%
% 2015-08-07 ~ v01.05.05:
% ~ add CDF_varInfo.m
% ~ use CDF FILLVALs
% ~ clarify use of (edi_BdvE_recnum:edi_xref2_BdvE) :: (1:many)
% ~ scale drift velocity from SI to km/s, IAW
%   https://lasp.colorado.edu/galaxy/display/mms/Units+of+Measure (2015-08-08)
%
% 2015-07-24 ~ v01.05.04:
% ~ add internal notes about global flags
% ~ replace references to 'target' with 'virtual source', or some variant
% ~ add myLibCDFConstants.m
%
% 2015-07-22 ~ v01.05.03:
% ~ add plot of weights versus intersection angles, with quality setpoints
% ~ add internal notes
% ~ change sinx_wt_Q_xovr definition; change plot to degrees for clarity
% ~ add beam_intercept_angle_test.m to project to validate calculation
% ~ add abs() to intersect angle - weight calcs (-90..90� ~> 0..90�)
%
% 2015-07-20 ~ v01.05.02:
% ~ add first estimate of beam intersection quality - parallelism
%
% ~ v01.05.01: IA < macroBeamCheckAngle = atan(tand(5)) ~> NaN
%   There are subtle (small, insignificant) differences in a few S*;
%   not worth the code expense.
% ~ v01.05.00: all beams included. sin^4 weighting
% ~ Use 'mms2_edi_slow_ql_efield_20150509_v0.1.3.cdf'
% ~ add intercept count to file name
% ~ add batch file to send project to Git dir
% ~ changed weight from sin^2 to sin^4
%
% 2015-06-30 ~ v01.04.00:
% ~ made plot figure persistent; removed figure close () in plot routine
% ~ finally got the plots to update & save in the background, but not sure that
%   there is an actual performance improvement.
% ~ add switch to choose among rotation matrix methods
% ~ add switch to show unit vectors
% ~ different symbols and lines for GDU1, 2
% ~ added legend to beam plots
%
% 2015-06-23 ~ v01.03.00:
% ~ included global vars to control program flow
% ~ made plots invisble during plotting if pause_each_plot = false
% ~ update trace disps to command window and plot title content
%
% 2015-06-22 ~ v01.02.00:
% ~ use CDF files for EDI, B, and EDP data
% ~ apply plot style during save
% ~ include obsID on plots
% ~ include local copies of myLibAppConstants, myLibScienceConstants
% ~ change references to Cluster data
%
% 2015-06-19 ~ v01.01.00:
% ~ update plot: sample datetime, selected confidence, legend
% ~ implement weighted mean: sin^2
% ~ test with 'mms2_edi_slow_ql_efield_20150509_v0.0.1.cdf'
%
% 2015-06-08 ~ v01.00.00:
% ~ draft, proof of concept
%

clc             % clear the command window
clear all
close all       % close all figures

format compact
format short g    % +, bank, hex, long, rat, short, short g, short eng
myLibAppConstants % custom colors; set default axis colors
myLibScienceConstants

global dotVersion;        dotVersion        = 'v1.07.00';
global pause_each_plot;   pause_each_plot   = false; % Plots are not visible if they are not paused
global plot_B_frunc;      plot_B_frunc      = false; % B fractional uncertainty
global plot_beams;        plot_beams        = false;
global plot_d_frunc;      plot_b_frunc      = false; % drift step fractional uncertainty
global plot_dots;         plot_dots         = false;
global plot_E_frunc;      plot_E_frunc      = false; % E fractional uncertainty
global plot_edi_E_sdcs;   plot_edi_E_sdcs   = false; % deprecated
global plot_summary;      plot_summary      = false;
global plot_sinx_weights; plot_sinx_weights = false;
global rotation_method;   rotation_method   = 2;     % 1: r_matrix, 2: Bx = By x Bz, 3: Euler x, y
global save_beam_plots;   save_beam_plots   = false;
sinx_wt_Q_xovr_angles                       = [ 8.0 30 ];
global sinx_wt_Q_xovr;    sinx_wt_Q_xovr    = sind (sinx_wt_Q_xovr_angles).^4.0; % breakpoints for quality ranges for sin^x weighting
global show_unit_vectors; show_unit_vectors = false;
global use_v10502;        use_v10502        = false; % Remove parallel beams

EDI_presentation_beam_plot_style.Resolution = '600';

UseFileOpenGUI = false; % true false
if ~UseFileOpenGUI
	mms_ql_dataPath = 'D:\MMS\MATLAB\MEEdrift\mms_edi_cdf';

% 	mms_ql_EDI_BdvE_dataFile = 'mms1_edi_srvy_sl_efield_20150819_v0.2.3.cdf';
% 	mms_ql_EDI_BdvE_dataFile = 'mms1_edi_srvy_sl_efield_20150820_v0.2.3.cdf';
% 	mms_ql_EDI_BdvE_dataFile = 'mms1_edi_srvy_sl_efield_20150821_v0.2.3.cdf';
% 	mms_ql_EDI_BdvE_dataFile = 'mms1_edi_srvy_sl_efield_20150822_v0.2.3.cdf';
% 	mms_ql_EDI_BdvE_dataFile = 'mms2_edi_srvy_sl_efield_20150819_v0.2.3.cdf';
% 	mms_ql_EDI_BdvE_dataFile = 'mms2_edi_srvy_sl_efield_20150820_v0.2.3.cdf';
% 	mms_ql_EDI_BdvE_dataFile = 'mms2_edi_srvy_sl_efield_20150821_v0.2.3.cdf';
% 	mms_ql_EDI_BdvE_dataFile = 'mms2_edi_srvy_sl_efield_20150822_v0.2.3.cdf';
	mms_ql_EDI_BdvE_dataFile = 'mms3_edi_srvy_sl_efield_20150819_v0.2.3.cdf';
% 	mms_ql_EDI_BdvE_dataFile = 'mms3_edi_srvy_sl_efield_20150820_v0.2.3.cdf';
% 	mms_ql_EDI_BdvE_dataFile = 'mms3_edi_srvy_sl_efield_20150821_v0.2.3.cdf';
% 	mms_ql_EDI_BdvE_dataFile = 'mms3_edi_srvy_sl_efield_20150822_v0.2.3.cdf';
% 	mms_ql_EDI_BdvE_dataFile = 'mms2_edi_srvy_sl_efield_20150820_v0.2.3.cdf';
% 	mms_ql_EDI_BdvE_dataFile = 'mms2_edi_srvy_sl_efield_20150821_v0.2.3.cdf';
% 	mms_ql_EDI_BdvE_dataFile = 'mms2_edi_srvy_sl_efield_20150822_v0.2.3.cdf';

	mms_ql_EDI_BdvE_data = [mms_ql_dataPath cFileSep mms_ql_EDI_BdvE_dataFile];

	mms_ql_dataPath = 'D:\MMS\MATLAB\MEEdrift\mms_edp_cdf';

% 	mms_ql_EDP_dataFile = 'mms1_edp_slow_ql_dce_20150819000000_v0.2.0.cdf';
% 	mms_ql_EDP_dataFile = 'mms1_edp_slow_ql_dce_20150820000000_v0.2.0.cdf';
% 	mms_ql_EDP_dataFile = 'mms1_edp_slow_ql_dce_20150821000000_v0.2.0.cdf';
% 	mms_ql_EDP_dataFile = 'mms1_edp_slow_ql_dce_20150822000000_v0.2.0.cdf';
% 	mms_ql_EDP_dataFile = 'mms2_edp_slow_ql_dce_20150819000000_v0.2.0.cdf';
% 	mms_ql_EDP_dataFile = 'mms2_edp_slow_ql_dce_20150820000000_v0.2.0.cdf';
% 	mms_ql_EDP_dataFile = 'mms2_edp_slow_ql_dce_20150821000000_v0.2.0.cdf';
% 	mms_ql_EDP_dataFile = 'mms2_edp_slow_ql_dce_20150822000000_v0.2.0.cdf';
	mms_ql_EDP_dataFile = 'mms3_edp_slow_ql_dce_20150819000000_v0.2.0.cdf';
% 	mms_ql_EDP_dataFile = 'mms3_edp_slow_ql_dce_20150820000000_v0.2.0.cdf';
% 	mms_ql_EDP_dataFile = 'mms3_edp_slow_ql_dce_20150821000000_v0.2.0.cdf';
% 	mms_ql_EDP_dataFile = 'mms3_edp_slow_ql_dce_20150822000000_v0.2.0.cdf';
% 	mms_ql_EDP_dataFile = 'mms4_edp_slow_ql_dce_20150820000000_v0.2.0.cdf';
% 	mms_ql_EDP_dataFile = 'mms4_edp_slow_ql_dce_20150821000000_v0.2.0.cdf';
% 	mms_ql_EDP_dataFile = 'mms4_edp_slow_ql_dce_20150822000000_v0.2.0.cdf';

	mms_ql_EDP_data = [mms_ql_dataPath cFileSep mms_ql_EDP_dataFile];
end

edi_drift_step_read_ql__EDI__B__EDP_data

nB_recs       = length (BdvE_dn);
d_bpp         = zeros (3, nB_recs);
d_SD_bpp  = zeros (3, nB_recs);
d_CI_bpp  = zeros (3, nB_recs);
d_sdcs        = zeros (3, nB_recs);
d_SD_sdcs = zeros (3, nB_recs);
d_CI_sdcs = zeros (3, nB_recs);
d_quality     = zeros (1, nB_recs, 'uint8');
v_sdcs  = zeros (3, nB_recs);
E_sdcs  = zeros (3, nB_recs);
E_sdcs_unc = zeros (3, nB_recs);

E_sdcs (:,:) = NaN;

% for each B interval, there are EDI records indexed to B:EDI :: 1:many
% Here we look at each B record, find the corresponding EDI records, and filter.

disp 'looping over edi_B_sdcs'
for B_recnum = 1: size (edi_B_sdcs, 2)
% 	sprintf (' processing B_recnum %d', B_recnum)
	% specific date-time request to compare w Matt's plot
	% MMS2 2015-05-09 16:08:35 - 16:08:40
	% 	tt2000_20150509_160835 = spdfdatenumtott2000(datenum('2015-05-09 16:08:34', 'yyyy-mm-dd HH:MM:ss'))
	% 	tt2000_20150509_160840 = spdfdatenumtott2000(datenum('2015-05-09 16:08:39', 'yyyy-mm-dd HH:MM:ss'))
	% 	B_recnum = find ( (edi_BdvE_t2k >= tt2000_20150509_160835) & (edi_BdvE_t2k < tt2000_20150509_160840) )

	BdvE_recnum  = edi_BdvE_recnum (B_recnum)'; % should be ONLY >1<
	iB_xref2_edi = find (edi_xref2_BdvE == BdvE_recnum); % 0 or more

	% keyboard
	if (length (iB_xref2_edi) > 1) % More than 1 beam
		B_sdcs    = edi_B_sdcs (:, B_recnum);
		B_SD_sdcs = edi_B_std_sdcs (:, B_recnum);
		B_t2k     = edi_BdvE_t2k (B_recnum);

		gd_virtual_sdcs = edi_gd_virtual_sdcs (:, iB_xref2_edi);
		gd_fv_sdcs      = edi_gd_fv_sdcs      (:, iB_xref2_edi);
		gd_ID           = edi_gd_ID           (:, iB_xref2_edi);

		[ ...
			d_bpp(:, B_recnum), ...
			d_SD_bpp(:, B_recnum), ...
			d_CI_bpp(:, B_recnum), ...
			d_sdcs(:, B_recnum), ...
			d_SD_sdcs(:, B_recnum), ...
			d_CI_sdcs(:, B_recnum), ...
			d_quality(:, B_recnum), ...
			v_sdcs(:, B_recnum), ...
			E_sdcs(:, B_recnum), ...
			E_sdcs_unc(:, B_recnum) ...
	  ] = edi_drift_step ( ...
			obsID, ...
			B_t2k, ...
			B_sdcs, ...
			B_SD_sdcs, ...
			gd_virtual_sdcs, ...
			gd_fv_sdcs, ...
			gd_ID ...
		);

	else
		E_sdcs (:, B_recnum) = [ NaN; NaN; NaN ];
	end

end % for B_recnum = 1: ...

% these next plots are debugging plots, and are meant to be used manually
% comment them out to run in usual mode
% keyboard

% plot (BdvE_dn, abs (edi_B_std_sdcs ./ edi_B_sdcs));
% title ('EDI B fractional uncertainty')

% plot (BdvE_dn, abs (d_SD_bpp ./ d_bpp));
% title ('edi\_drift\_step d BPP stddev fractional uncertainty')

% plot (BdvE_dn, abs (d_CI_bpp ./ d_bpp));
% title ('edi\_drift\_step d BPP uncertainty from mean stddev fractional uncertainty')

% plot (BdvE_dn, abs (d_SD_sdcs ./ d_sdcs));
% title ('edi\_drift\_step d SDCS stddev fractional uncertainty')

% plot (BdvE_dn, abs (d_CI_sdcs ./ d_sdcs));
% title ('edi\_drift\_step d SDCS uncertainty fractional uncertainty')

% plot (BdvE_dn, abs(E_sdcs_unc ./ E_sdcs));
% title ('E\_sdcs\_unc fractional uncertainty')

% plot (BdvE_dn, E_sdcs, BdvE_dn, edi_E_sdcs);
% legend ('E sdcs_x', 'E sdcs_y', 'E sdcs_z', 'EDI E_x', 'EDI E_y', 'EDI E_z');
% title ('E SDCS and EDI E')

%{
xlabel (sprintf ( '%s', datestr((BdvE_dn(1)), 'yyyy-mm-dd') ))
plotDateMin = double (min(BdvE_dn));
plotDateMax = double (max(BdvE_dn));
% plotDateMin = datenum ('2015-08-19 09:30:00'); % datestr (plotDateMin, 'yyyy-mm-dd HH:MM:SS');
plotDateMax = datenum ('2015-08-19 09:40:00');
xlim ( [ plotDateMin plotDateMax ] )
set (gca, 'XTick', [ plotDateMin: datenum_10min: plotDateMax])
datetick ('x', 'HH:MM', 'keeplimits', 'keepticks')
%}

save edi_drift_step_20150819_105000_105500_all_vars.mat
% load edi_drift_step_20150819_105000_105500_all_vars.mat
% End of manual debugging plots - comment them out to run in usual mode

if plot_summary
	set (0, 'DefaultFigureVisible', 'on')

	[ hAxis hEDP_dce_xyz_sdcs hEDI_B ] = plotyy ( ...
		edp_dn, edp_dce_xyz_sdcs(:, 1:2), ...
		BdvE_dn, edi_B_sdcs (1:3, :), @plot, @plot );

	% retained for full-width data plot
	plotDateMin = double (min (min(edp_dn), min(BdvE_dn)));
	plotDateMax = double (max (max(edp_dn), max(BdvE_dn)));

	% 	plotDateMin = datenum ('2015-05-09 15:40:00'); % datestr (plotDateMin, 'yyyy-mm-dd HH:MM:SS');
	% 	plotDateMax = datenum ('2015-05-09 16:40:00');

	xlim ( [ plotDateMin plotDateMax ] )
	xlabel (sprintf ( '%s', datestr( spdftt2000todatenum (edp_t2k(1)), 'yyyy-mm-dd') ))
	set (gca, 'XTick', [ plotDateMin: datenum_1hr: plotDateMax])
	datetick ('x', 'HH:MM', 'keeplimits', 'keepticks')
	set (hAxis(2), 'XTick', [])

% 	ylim ([ -4 4 ])
	ylabel (hAxis (1),'mV-m^-^1')
	ylabel (hAxis (2),'nT')

	grid on
	hold on

	EDI_E_B_plotColors = [ ...
		myDarkTeal; myGrassGreen; ...
		myDarkBlue; myDarkGreen; ...
		myLightGrey6; ...
		MMS_plotColorx; MMS_plotColory; MMS_plotColorz ];

	set (hEDP_dce_xyz_sdcs(1), 'Color', myDarkTeal, ...
		'LineStyle', 'none', 'Marker', 'o', ...
		'MarkerFaceColor', myDarkTeal, 'MarkerEdgeColor', myDarkTeal, 'MarkerSize', 2.0);
	set (hEDP_dce_xyz_sdcs(2), 'Color', myGrassGreen, ...
		'LineStyle', 'none', 'Marker', 'o', ...
		'MarkerFaceColor', myGrassGreen, 'MarkerEdgeColor', myGrassGreen, 'MarkerSize', 2.0);

	set (hEDI_B(1), 'Color', MMS_plotColorx, ...
		'LineStyle', 'none', 'Marker', 'o', ...
		'MarkerFaceColor', MMS_plotColorx, 'MarkerEdgeColor', MMS_plotColorx, 'MarkerSize', 1.0);
	set (hEDI_B(2), 'Color', MMS_plotColory, ...
		'LineStyle', 'none', 'Marker', 'o', ...
		'MarkerFaceColor', MMS_plotColory, 'MarkerEdgeColor', MMS_plotColory, 'MarkerSize', 1.0);
	set (hEDI_B(3), 'Color', MMS_plotColorz, ...
		'LineStyle', 'none', 'Marker', 'o', ...
		'MarkerFaceColor', MMS_plotColorz, 'MarkerEdgeColor', MMS_plotColorz, 'MarkerSize', 1.0);

	hE_drift_sdcsx = plot (BdvE_dn, E_sdcs (1,:), ...
		'LineStyle', 'none', 'Marker', 'o', ...
		'MarkerFaceColor', myDarkBlue, 'MarkerEdgeColor', myDarkBlue, 'MarkerSize', 3.0);

	hE_drift_sdcsy = plot (BdvE_dn, E_sdcs (2,:), ...
		'LineStyle', 'none', 'Marker', 'o', ...
		'MarkerFaceColor', myDarkGreen, 'MarkerEdgeColor', myDarkGreen, 'MarkerSize', 3.0);

% 	hQuality_plot = bar (BdvE_dn, d_quality, 0.1, 'w', 'EdgeColor', 'white'); % myLightGrey6

	hEDI_E_B_plot (1:2) = hEDP_dce_xyz_sdcs;
	hEDI_E_B_plot (3:4) = [ hE_drift_sdcsx, hE_drift_sdcsy ];
% 	hEDI_E_B_plot (5)   = hQuality_plot;
	hEDI_E_B_plot (6:8) = hEDI_B;

	title ([' EDI drift step ', dotVersion, ': MMS', obsID, ' SDP 2D E-fields, EDI E-fields and B avg @ 5 s intervals, DMPA'])

	hLegend = legend (hEDI_E_B_plot, 'SDP E_x', 'SDP E_y', 'EDI E_x', 'EDI E_y', 'Quality', 'B_x', 'B_y', 'B_z');
	hText = findobj (hLegend, 'type', 'text');
	nhText = size (hText,1);
	for i = 1: nhText
		set (hText (nhText-i+1), 'color', EDI_E_B_plotColors (i,:)); % Believe it or not, hText is in reverse order!!
	end

	hold off

	figure
	hist (double(d_quality), [0:3])
end % if plot_summary

if plot_sinx_weights
	figure;
	hAxes = axes;
	intersect_angles = 0.0: 0.1: 90.0;
	plot (intersect_angles, sind(intersect_angles).^4.0)
	xlim ([0.0 90.0])
	hold on
	title 'Sin^4 weights versus intersection angles, with quality setpoints'
	line ([sinx_wt_Q_xovr_angles(1) sinx_wt_Q_xovr_angles(1)], get (hAxes, 'YLim'), 'Color', [1.0 0.0 0.0]);
	line ([sinx_wt_Q_xovr_angles(2) sinx_wt_Q_xovr_angles(2)], get (hAxes, 'YLim'), 'Color', [1.0 0.0 0.0]);
	text ( 3.0, 0.5, 'Q1')
	text (17.0, 0.5, 'Q2')
	text (60.0, 0.5, 'Q3')
	xlabel ('Intercept angle, degrees')
	ylabel ('Weight')
	hold off
end % if plot_sinx_weights
