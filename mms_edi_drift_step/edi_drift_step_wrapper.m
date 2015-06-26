% edi_drift_step using beam convergence, v0101 (show also on plots)

% v0101: update plot: sample datetime, selected confidence, legend
%        implement weighted mean: sin^2
%        tested with 'mms2_edi_slow_ql_efield_20150509_v0.0.1.cdf'
% v0100: draft, proof of concept

clc             % clear the command window
clf
clear variables
close all       % close all figures
format compact
format short g  % +, bank, hex, long, rat, short, short g, short eng
myLibAppConstants % custom colors; set default axis colors

edi_drift_step_read_ql__EDI__B__data % Get all the data from (ex.) mms2_edi_slow_ql_efield_20150509_v0.0.1.cdf

% for each B interval, there are EDI records indexed to B ~ many EDI:1 B
% Here we look at each B record, find the corresponding EDI records, and filter.
for igdxx = 1: length (edi_B_dmpa)
disp 'looping over edi_B_dmpa'
	B_dmpa   = edi_B_dmpa      (1:3,igdxx);
	if ~isnan (B_dmpa (1))

		% B field data
		B_tt2000 = edi_BdvE_tt2000 (    igdxx);

		% GDU data that corresponds to the B field data: position and firing vectors
		% Includes both GDUs
		iigd12_b_avgIntrp     = find (edi_gd_B_index == igdxx);

		gd_virtual_dmpa       = edi_gd_virtual_dmpa (:, iigd12_b_avgIntrp);
		gd_fv_dmpa            = edi_gd_fv_dmpa      (:, iigd12_b_avgIntrp);
		gd_ID                 = edi_gd_ID           (:, iigd12_b_avgIntrp);

		% keyboard
		if (length (gd_virtual_dmpa) > 2)
			E_dmpa (:, igdxx) = edi_drift_step ( ...
				B_tt2000, ...
				B_dmpa, ...
				gd_virtual_dmpa, ...
				gd_fv_dmpa, ...
				gd_ID );
		else
			E_dmpa (:, igdxx) = [ NaN; NaN; NaN ];
		end

	end
end
save 'E_dmpa_v0101.mat' E_dmpa

%{
edp_tt2000 = edp_ql.tt2000; % Epoch times
e_dsl      = edp_ql.e_dsl;  % Electric field in DSL coordinates (DSL ~= DMPA)

edp_epoch  = double (edp_tt2000) * 1e-9; % ~> sec
bAvg_epoch = double (b_avg.t_avg) * 1e-9;
[ datestr(spdftt2000todatenum(edp_tt2000(1)), 'yyyy-mm-dd HH:MM:ss'), ' ',...
  datestr(spdftt2000todatenum(edp_tt2000(end)), 'yyyy-mm-dd HH:MM:ss') ]
edp_datenum  = spdftt2000todatenum(edp_tt2000');
bAvg_datenum = spdftt2000todatenum(b_avg.t_avg');

[ hAxis hE_dsl hB_avg ] = plotyy ( ...
	edp_datenum, e_dsl(1:2,:)', ...
	bAvg_datenum, b_avg.b_avg (1:3,:)' );

datenumOneMin = 0.0006944444496185;
xlim ( [ datenum('2015-05-06 15:30:00') datenum('2015-05-06 15:35:00') ] )
 set (gca, 'XTick', [datenum('2015-05-06 15:30:00'): datenumOneMin: datenum('2015-05-06 15:35:00')])

datetick ('x', 'HH:MM', 'keeplimits', 'keepticks')

set (hAxis(2), 'XTick', [])

hold on
% set (hB_avg, 'Color', [ MMS_plotColorx, MMS_plotColory, MMS_plotColorz ]);
set (hB_avg(1), 'Color', MMS_plotColorx);
set (hB_avg(2), 'Color', MMS_plotColory);
set (hB_avg(3), 'Color', MMS_plotColorz);

plot (bAvg_datenum, E_dmpa (1,:), ...
	'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', MMS_plotColorx, 'MarkerEdgeColor', MMS_plotColorx, 'MarkerSize', 5.0);
% datetick ()

plot (bAvg_datenum, E_dmpa (2,:), ...
	'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', MMS_plotColory, 'MarkerEdgeColor', MMS_plotColory, 'MarkerSize', 5.0);
% datetick ()
title 'SDP 2D E-fields, EDI E-fields and B avg @ 5 s intervals, DMPA'
xlabel (sprintf ( '%s', datestr( spdftt2000todatenum (edp_tt2000(1)), 'yyyy-mm-dd') ))
ylabel (hAxis (1),'mV-m^-^1')
ylabel (hAxis (2),'nT')
legend ('E_x SDP', 'E_y SDP', 'E_x EDI', 'E_y EDI', 'B_x', 'B_y', 'B_z' );

hold off
%}

% b_avg_tt2000 = b_avg.t_avg; % MATLAB does not save structure.variables
% save 'mms4_edi_slow_l1a_efield_20150506_SDP_and_EDI_driftstep_E_field_2D.mat' edp_tt2000 e_dsl b_avg_tt2000 E_dmpa
% load 'mms4_edi_slow_l1a_efield_20150506_SDP_and_EDI_driftstep_E_field_2D.mat'

%{

Some of the variables in the mat file
% Restore the data
mat_file = 'mms4_edi_slow_l1a_efield_20150506_v0.1.0.mat';
load (mat_file)

% Average Quantities
tt2000_avg = b_avg.t_avg;        % time tags of 5-second avg B-field
dt_avg     = b_avg.dt_avg;       % delta plus of t_avg
b_average  = b_avg.b_avg;        % Avg B-field
b_std      = b_avg.b_std;        % Standard deviation of b_avg
b_gd12     = b_avg.b_gd12;       % B interpolated to Gun1 times
b_gd21     = b_avg.b_gd21;       % B interpolated to Gun2 times
inds_gd12  = b_avg.inds_gd12;    % Array of indices indicating to which b_avg value a firing vector maps.
inds_gd21  = b_avg.inds_gd21;    % Array of indices indicating to which b_avg value a firing vector maps.

% EDI data (the structure has more)
tt2000_gd12       = edi.epoch_gd12;              % Epoch times of gun1
tt2000_gd21       = edi.epoch_gd21;              % Epoch times of gun2
virtual_gun1_dmpa = edi.virtual_gun1_dmpa;       % Location of gun1 on virtual spacecraft
virtual_gun2_dmpa = edi.virtual_gun2_dmpa;       % Location of gun2 on virtual spacecraft
fv_gd12_dmpa      = edi.fv_gd12_dmpa;            % Firing vectors from gun1
fv_gd21_dmpa      = edi.fv_gd21_dmpa;            % Firing vectors from gun2

% FG Data
fg_tt2000 = fg_ql.tt2000;       % Epoch times
b_dmpa    = fg_ql.b_dmpa;       % Magnetic field

% EDP Data
edp_tt2000 = edp_ql.tt2000;     % Epoch times
e_dsl      = edp_ql.e_dsl;      % Electric field in DSL coordinates (DSL ~= DMPA)

%}