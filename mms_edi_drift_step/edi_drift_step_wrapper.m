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
% 'mms2_edi_slow_ql_efield_20150509_v0.0.1.cdf'
% 'mms2_edp_comm_ql_dce2d_20150509000000_v0.0.0.cdf'
%
% Output:
% Optionally produces a series of drift step plots.
% Optionally saves drift step plots.
% Optionally produces a summary plot of B, SDP, and drift step E-field.
%
% MATLAB release(s) MATLAB 8.3.0.532 (R2014a)
% Required Products None
%
% History:
% 2015-06-xx ~ v0104:
%  ~ made plot figure persistent; removed figure close () in plot routine
%  ~ finally got the plots to update & save in the background, but not sure that
%    there is an actual performance improvement.
%  ~ added switch to choose among rotation matrix methods
%  ~ added switch to show unit vectors
%  ~ different symbols and lines for GDU1, 2
%  ~ added legend to beam plots
% 2015-06-23 ~ v0103:
%  ~ included global vars to control program flow
%  ~ made plots invisble during plotting if pause_each_plot = false
%  ~ updated trace disps to command window and plot title content
% 2015-06-22 ~ v0102:
%  ~ use CDF files for EDI, B, and EDP data
%  ~ apply plot style during save
%  ~ include obsID on plots
%  ~ include local copies of myLibAppConstants, myLibScienceConstants
%  ~ changed references to Cluster data
% 2015-06-19 ~ v0101:
%  ~ update plot: sample datetime, selected confidence, legend
%  ~ implement weighted mean: sin^2
%  ~ tested with 'mms2_edi_slow_ql_efield_20150509_v0.0.1.cdf'
% 2015-06-08 ~ v0100:
%  ~ draft, proof of concept
%

clc             % clear the command window
clear variables
close all       % close all figures
format compact
format short g  % +, bank, hex, long, rat, short, short g, short eng
myLibAppConstants % custom colors; set default axis colors

global dotVersion;        dotVersion        = 'v1.04';
global pause_each_plot;   pause_each_plot   = true;   % Plots are not visible if they are not paused
global plot_beams;        plot_beams        = true;
global plot_dots;         plot_dots         = false;
global plot_summary;      plot_summary      = true;
global rotation_method;   rotation_method   = 2;      % 1: r_matrix, 2: Bx = By x Bz, 3: Euler x, y
global save_beam_plots;   save_beam_plots   = true;
global show_unit_vectors; show_unit_vectors = true;

edi_drift_step_read_ql__EDI__B__EDP_data

nB_recs       = length (B_datenum);
driftStep     = zeros (3, nB_recs);
dUncertainty  = zeros (3, nB_recs);
driftVelocity = zeros (3, nB_recs);
drift_E       = zeros (3, nB_recs);
drift_E (:,:) = NaN;

% for each B interval, there are EDI records indexed to B ~ many EDI:1 B
% Here we look at each B record, find the corresponding EDI records, and filter.
disp 'looping over edi_B_dmpa'
for B_recnum = 1: length (edi_B_dmpa)

	% specific date-time request to compare w Matt's plot
	% MMS2 2015-05-09 16:08:35 - 16:08:40
	% 	tt2000_20150509_160835 = spdfdatenumtott2000(datenum('2015-05-09 16:08:34', 'yyyy-mm-dd HH:MM:ss'))
	% 	tt2000_20150509_160840 = spdfdatenumtott2000(datenum('2015-05-09 16:08:39', 'yyyy-mm-dd HH:MM:ss'))
	% 	B_recnum = find ( (edi_BdvE_tt2000 >= tt2000_20150509_160835) & (edi_BdvE_tt2000 < tt2000_20150509_160840) )

	BdvE_recnum  = edi_BdvE_recnum (B_recnum)'; % should be ONLY >1<
	iB_to_edi = find (edi_gd_B_xref == BdvE_recnum);

	B_dmpa   = edi_B_dmpa (1:3,B_recnum);
	if ~isnan (B_dmpa (1))

		% B field data
		B_tt2000 = edi_BdvE_tt2000 (B_recnum);

		% keyboard
		if (length (iB_to_edi) > 2) % More than 2 beams

			gd_virtual_dmpa = edi_gd_virtual_dmpa (:, iB_to_edi);
			gd_fv_dmpa      = edi_gd_fv_dmpa      (:, iB_to_edi);
			gd_ID           = edi_gd_ID           (:, iB_to_edi);

			[ driftStep(:, B_recnum), ...
			  dUncertainty(:, B_recnum), ...
			  driftVelocity(:, B_recnum), ...
			  drift_E(:, B_recnum) ] = edi_drift_step ( ...
				B_tt2000, ...
				B_dmpa, ...
				gd_virtual_dmpa, ...
				gd_fv_dmpa, ...
				obsID, ...
				gd_ID );
		else
			drift_E (:, B_recnum) = [ NaN; NaN; NaN ];
		end

	end
end % for B_recnum = 1: ...

if plot_summary
	set (0, 'DefaultFigureVisible', 'on')

	[ hAxis hEDP_dce_xyz_dsl hEDI_B ] = plotyy ( ...
		edp_datenum, edp_dce_xyz_dsl(:, 1:2), ...
		B_datenum, edi_B_dmpa (1:3, :)', @plot, @plot );

	plotDateMin = double (min (min(edp_datenum), min(B_datenum)));
	plotDateMax = double (max (max(edp_datenum), max(B_datenum)));

	datenumOneMin = 0.0006944444496185;
	datenumOneHr  = 60.0 * datenumOneMin;
	xlim ( [ plotDateMin plotDateMax ] )
	set (gca, 'XTick', [ plotDateMin: datenumOneHr: plotDateMax]) % debug?

	datetick ('x', 'HH:MM', 'keeplimits', 'keepticks') % debug?

	set (hAxis(2), 'XTick', [])

	hold on
	set (hEDI_B(1), 'Color', MMS_plotColorx);
	set (hEDI_B(2), 'Color', MMS_plotColory);
	set (hEDI_B(3), 'Color', MMS_plotColorz);

	hDrift_Ex = plot (B_datenum, drift_E (1,:), ...
		'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', MMS_plotColorx, 'MarkerSize', 5.0);

	hDrift_Ey = plot (B_datenum, drift_E (2,:), ...
		'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', myDarkGreen, 'MarkerEdgeColor', MMS_plotColory, 'MarkerSize', 5.0);
	% datetick ()
	title 'SDP 2D E-fields, EDI E-fields and B avg @ 5 s intervals, DMPA'
	xlabel (sprintf ( '%s', datestr( spdftt2000todatenum (edp_tt2000(1)), 'yyyy-mm-dd') ))
	ylabel (hAxis (1),'mV-m^-^1')
	ylabel (hAxis (2),'nT')
	legend ('SDP E_x', 'SDP E_y', 'EDI E_x', 'EDI E_y', 'B_x', 'B_y', 'B_z' );

	hold off
end
