% MEEdrift
MEEdrift_header % Please see the header for important application notes and programming conventions

if ~exist ('Reentrant')
	Reentrant = false;
	strReentrant = 'Reentrancy is OFF';
end
if ~Reentrant
	disp 'Initializing...'
	clc             % clear the command window
	clear variables
	close all       % close all figures
	Reentrant = false; % must be redeclared because it was just cleared
	strReentrant = 'Reentrancy is OFF';
end

if ~Reentrant % If Reentrant = true, don't run startup dialogs and don't initialize global vars
	global Project_;       Project_       = 'MMS_EDI_EDP_driftstep';
	global dotVersion;     dotVersion     = 'v1.00.00';
	global beamsInWindow;  beamsInWindow  = 8; % n, NOTE !!!! ± beamsInWindow before|after center beam time
	global beamWindowSecs; beamWindowSecs = 4; % (s), NOTE !!!! ± beamWindowSecs before|after center beam time
	beamWindowDays = beamWindowSecs / 86400.0; % Convert from seconds to days; to compare with datenums
	global PlotHoldOff;    PlotHoldOff    = true;  % Set to false to plot multiple beam data records on the same plot; true gets but 1            
	global PlotSubplots;   PlotSubplots   = false; % true=plot the data as it is read in; only works for read routines with embedded plot commands
	global PlotView;       PlotView       = '2D';  % or 3D; the plot is always 3D, but this switches the initial plot display perspective         
	global ScrollData;     ScrollData     = true;
	global UseFileOpenGUI; UseFileOpenGUI = true;
	global Use_OCS_plot;   Use_OCS_plot   = false;

	RichTest = false; % true false

	myLibAppConstants
	% Constants used to calculate beam geometry
	myLibScienceConstants

	hEDP_plot  = 1;
	hDMPA_plot = 2;
	hBPP_plot  = 3;
	hFFT_plot  = 4;

	% [left, bottom, width, height]
	hEDP_plot = figure ('Position', [ 0.0*ScreenWidth 0.55*ScreenHeight 0.9*ScreenWidth 0.375*ScreenHeight ]);
	set (hEDP_plot, 'WindowStyle', 'normal')
	set (hEDP_plot, 'DockControls', 'off')
	set (gcf, 'name', 'MMS EDP Ex Ey (mV)', 'NumberTitle', 'off');

	hDMPA_plot = figure ('Position', [ 0.0*ScreenWidth 0.0*ScreenHeight 0.7*ScreenHeight 0.6*ScreenHeight ]);
	set (hDMPA_plot, 'WindowStyle', 'normal')
	set (hDMPA_plot, 'DockControls', 'off')
	set (gcf, 'name', 'Spacecraft and beams in DMPA', 'NumberTitle', 'off');

	hBPP_plot = figure ('Position', [ 0.5*ScreenWidth 0.0*ScreenHeight 0.7*ScreenHeight 0.6*ScreenHeight ]);
	set (hBPP_plot, 'WindowStyle', 'normal')
	set (hBPP_plot, 'DockControls', 'off')
	set (gcf, 'name', 'Spacecraft and beams in BPP', 'NumberTitle', 'off');

	% MMS Phase 1 orbits (19 000 to 119 000 kilometres from the planet)
	% MMS speed @ perigee @ 19000 km ~= 3964 m/s    @ apogee @ 119000 km ~= 1784 m/s
	% Distance spacecraft moves in gyroPeriod of 0.00023835 s @ 3000 m/s ~= 0.715 m
	% MMS spacecraft mass at launch 1,200 kilograms

	% The flight software on MMS adjusts the beam firing direction within BPP
	% by 0.9? degrees every 3.9? milliseconds. The numbers are approximate, but good enough
	% for the purpose of adjusting the beam angular uncertainty for multi-runners.
	% @ 0.05 Hz -> 360°/20s -> 18°/s -> ??0.3515°/3.906 ms
	% - actual stepping size: ??arctan (1/64) =~ 0.895 degrees in 0.015625 s
	% - processing loop period: ??1/256 seconds =~ 3.906 ms
	mmsRPM           = 3.0; % nominal, 1 spin ~20 s. 0.05 Hz
	mmsSpinPeriod    = 20.0; % seconds, nominal
	EDP_sampleFreq   = 128.0; % samples per second
	EDP_samplePeriod = 1.0 / 128.0; % seconds per sample
	GDUz_loc = 0.0;

	% See 'MMS Hardware and Data Processing notes.txt'
	EDI1gunLoc = [ -1.45598,  1.11837, 0.0 ]; % EDI2 gun atan2(-1.11837, 1.45598)*180/pi ~> -37.52865°
	EDI1detLoc = [ -1.35885,  1.03395, 0.0 ]; % EDI2 det atan2(-1.03395, 1.35885)*180/pi ~> -37.26753°
	EDI2gunLoc = [  1.45598, -1.11837, 0.0 ]; % EDI1detLoc:EDI1gunLoc angle = atan2(1.11837-1.03395, -1.45598+1.35885)*180/pi ~> -40.995°
	EDI2detLoc = [  1.35885, -1.03395, 0.0 ]; % norm(EDI1gunLoc-EDI1detLoc,2) = 0.128689
	mmsEDI_VirtualRadius = norm (EDI1gunLoc-EDI2detLoc, 2);

	instrPlaneDMPA    = zeros (3, 361, 'double');
	for theta = 0:1:360
		instrPlaneDMPA (:, theta+1) = mmsEDI_VirtualRadius * [ cosd(theta); sind(theta); 0.0 ];
	end

	BeamCoordsX  = -8: 0.01: 8;
	Beam_z       = zeros (1, size (BeamCoordsX, 2), 'double');

	% default 1; ii = index into index, a doubly-dereferenced pointer
	% iiSorted_beam_t2k is used to scroll through iSorted_beam_t2k in linear
	% order, and iSorted_beam_t2k is an index of time-sorted pointers into *beam_t2k
	iiSorted_beam_t2k = uint32 (1);
	FFTplotted = false;
end

disp (['Starting MEEdrift ' dotVersion])

ButtonReply = '';
if ~Reentrant % If Reentrant = true, don't run startup dialogs and don't initialize global vars
	DialogPos = [ 600 500 ];
	[ButtonReply, DialogPos] = MEEdrift_questDlg_uitools (DialogPos, ...
		'Load Saved or Process New Data', 'Data Options', 'Load Saved', 'Process New', 'Process New'); % limit 3 buttons
	switch ButtonReply
		case 'Load Saved'
			% ==================================== Choose your time period and spacecraft
			[eventFile, eventPath] = uigetfile ('MEEdrift*.mat', 'Select date, spacecraft, time to analyze');
			if isequal (eventFile,  0) % then no valid file selected
				msgbox ('No valid MEEdrift event file selected.');
				ValidDataLoaded = false;
			else
				load ([eventPath, eventFile])
			end
			LoadEvent = true;
			SaveEvent = false;
		case 'Process New'
			LoadEvent = false;
			SaveEvent = false;
	end

	if LoadEvent
	else
		if RichTest
			beamsInWindow = 8;
			UseFileOpenGUI = false;

			mms_ql_dataPath = 'D:\MMS\MATLAB\MEEdrift';

			mms_ql__EDI__BdvE__dataFile = 'mms2_edi_slow_ql_efield_20150509_v0.1.4.cdf';
			mms_ql__EDI__BdvE__data = [mms_ql_dataPath cFileSep mms_ql__EDI__BdvE__dataFile];

			mms_ql__EDP__dataFile = 'mms2_edp_comm_ql_dce2d_20150509120000_v0.1.0.cdf';
			mms_ql__EDP_data = [mms_ql_dataPath cFileSep mms_ql__EDP__dataFile];

			Event = mms_ql__EDI__BdvE__dataFile (1:24); % An event label used to indentify this analysis
			Event = ['MEEdrift_', dotVersion, '_', Event];
		else
			ValidDataLoaded = false;
			UseFileOpenGUI = true;
		end

		MEEdrift_read_ql__EDI__B__EDP_data
	end
else
	ValidDataLoaded = true; % assume that user knows data is good
end % if ~Reentrant

ButtonReply  = '';
DSI_view3D = true;
Selection = 2;
strSelection = '> Beam >';
plotBeamDots = false;
strPlotBeamDots = 'Plot Beam Dots is OFF';

if ScrollData
	while ValidDataLoaded && ~strcmp (strSelection, 'Quit')
		if ~isempty (strfind ([ ...
			'< Beam <', ...
			'> Beam >', ...
			'Skip via EDP', ...
			'Plot Beam Dots is OFF', ...
			'Plot Beam Dots is ON' ], ...
			strSelection))
			MEEdrift_plot_beams
		end

		Default_UIControlFontSize = get (0, 'DefaultUIControlFontSize');
		set (0, 'DefaultUIControlFontSize', 12)
		ListCaptions = { ...
			'< Beam <', ...
			'> Beam >', ...
			strReentrant, ...
			'Skip via EDP', ...
			'Save Plots', ...
			'Save Snapshot', ...
			'Toggle DMPA <> BPP', ...
			'FFT plots of EFW x,y', ...
			strPlotBeamDots};

		[Selection, OK] = listdlg ( ...
			'Name', '', ...
			'PromptString', 'Select One Operation', ...
			'SelectionMode', 'single', ...
			'CancelString', 'Quit', ...
			'InitialValue', [Selection], ...
			'ListSize', [ 200 200 ], ...
			'ListString', ListCaptions);

		set (0, 'DefaultUIControlFontSize', Default_UIControlFontSize);
		if OK == 1
			strSelection = ListCaptions {Selection};
			switch strSelection
			  case '< Beam <'
					iiSorted_beam_t2k = max (1, iiSorted_beam_t2k-1);

				case '> Beam >'
					iiSorted_beam_t2k = min (nBeams, iiSorted_beam_t2k+1);

				case 'Reentrancy is OFF'
					Reentrant = true;
					strReentrant = 'Reentrancy is ON';

				case 'Reentrancy is ON'
					Reentrant = false;
					strReentrant = 'Reentrancy is OFF';

				case 'Skip via EDP'
					figure (hEDP_plot);
					[edp_timeIndex, mV] = ginput (1); % time axis is in datenum, E-field in mV
					last_iiSorted_beam_t2k = iiSorted_beam_t2k;
					if ~isempty (edp_timeIndex)
						% Go to the nearest EDI record that is in the neighborhood of this EDP record
						% Note the switch from EDP to EDI here
% disp ( sprintf ('Skip via EDP: edp_timeIndex: s%', datestr (edp_timeIndex, 'yyyy-mm-dd HH:MM:ss.fff') ) ) % V&V
						iedi_timeIndex = find (edi_gd_beam_dn > edp_timeIndex, 1);
% disp ( sprintf ('Skip via EDP: edi_gd_beam_dn %s', datestr (edi_gd_beam_dn (iedi_timeIndex), 'yyyy-mm-dd HH:MM:ss.fff') ) ) % V&V
						if ~isempty (iedi_timeIndex)
							iiSorted_beam_t2k = find (iSorted_beam_t2k == iedi_timeIndex, 1);
						else
							disp 'Warning!!! No beam data found for this EDP point!'
							iiSorted_beam_t2k = last_iiSorted_beam_t2k;
						end
					end

				% sample filename: 'MEEdrift_1.00.00_M2_20010305_060419.217@20150108.224134*'
				% mms2_edi_slow_ql_efield_20150509_v0.1.4.cdf
				case 'Save Plots'
					SavePlotFilename = [ ...
						'.' cFileSep 'MEEdrift_', dotVersion, ...
						'_M', obsID, mms_ql__EDP__dataFile(23:31), '_', ...
						strBeamTime(12:13), strBeamTime(15:16), strBeamTime(18:23), '_', ...
						datestr(now, '@yyyymmdd_HHMMSS')];
					saveas (hEDP_plot,  [ SavePlotFilename, 'a.png' ], 'png');
					if Use_OCS_plot
						saveas (hDMPA_plot, [ SavePlotFilename, 'b.png' ], 'png');
					end
					saveas (hBPP_plot,  [ SavePlotFilename, 'c.png' ], 'png');
					if FFTplotted
						saveas (hFFT_plot, [ SavePlotFilename, 'd.png' ], 'png');
					end
					disp (['Images saved as     "', SavePlotFilename, '[a|b|c|d].png', '".'])
					MEEdrift_writeCurrentPlotData

				case 'Save Snapshot'
				  disp ([ 'Event: ', Event ])
					if exist (Event, 'file')
						[status, attribs] = fileattrib ([Event,'.mat'])
						WriteMatFileOK = ~attribs.UserRead
					else
						WriteMatFileOK = true;
					end
					if WriteMatFileOK
						save (Event)
					else
						% 			warndlg (['Cannot save: ', Event, '.mat is Read Only.'], '!! Warning !!', 'modal') % does not work as modal r2014a
						% 			errordlg (['Cannot save: ', Event, '.mat is Read Only.'], '!! Warning !!', 'modal') % does not work as modal r2014a
						questdlg (['Cannot save: ', Event, '.mat is Read Only.'], 'Warning', 'OK', 'Cancel', 'OK'); % limit 3 buttons
					end

				case 'Toggle DMPA <> BPP'
					DSI_view3D = ~DSI_view3D;
					if DSI_view3D
						if Use_OCS_plot
							figure (hDMPA_plot);
							view ([ 115 20 ]); % Azimuth, Elevation in degrees, 3D view
						end
					else % 2D view of BPP
						% Azimuth, Elevation in degrees, 3D view, based on B DMPA vector
						figure (hBPP_plot);
						view ([   0 90 ]); % Azimuth, Elevation in degrees, std x-y plot
						% view ([ 42 58 ])
						%	view ([ atand(abs(centerBeamB_u(2)/centerBeamB_u(1))) acosd(centerBeamB_u(3))-90.0 ]);
					end

				case 'FFT plots of EFW x,y'
					MEEdrift_FFT_MMS_EFW_data

				case 'Plot Beam Dots is OFF'
					plotBeamDots = true;
					strPlotBeamDots = 'Plot Beam Dots is ON';

				case 'Plot Beam Dots is ON'
					plotBeamDots = false;
					strPlotBeamDots = 'Plot Beam Dots is OFF';

			end % switch strSelection
		else % if OK == 1
			strSelection = 'Quit';
		end
	end % while ValidDataLoaded && ~strcmp (strSelection, 'Quit')

	hold off
	% clear gyroFrequency GyroPeriod_SI BavgBPP d_PredExB;
else % if ScrollData; else batch process
	MEEdrift_batch_process
end

disp (['MEEdrift ' dotVersion ' ended'])

MEEdrift_references
