% edi_drift_step_read_ql__EDI__B__EDP_data

UseFileOpenGUI = true;
myLibCDFConstants

% ~~~~~~~~~~~~~~~~~~~
% mms2_edi_slow_ql_efield_20150509_v0.1.4.cdf'
% ~~~~~~~~~~~~~~~~~~~
if UseFileOpenGUI
	[mms_ql__EDI__BdvE__dataFile, mms_ql_dataPath] = uigetfile ('mms*.cdf', 'Select an MMS ql EDI_&_B CDF file');
	if isequal (mms_ql__EDI__BdvE__dataFile,  0) % then no valid file selected
		msgbox ('No valid MMS ql EDI_&_B data file selected.');
	else
		mms_ql__EDI__BdvE__data = [mms_ql_dataPath, mms_ql__EDI__BdvE__dataFile];
	end
end

obsID = mms_ql__EDI__BdvE__dataFile (4:4);
YYYY  = str2num (mms_ql__EDI__BdvE__dataFile (30:33));
MM    = str2num (mms_ql__EDI__BdvE__dataFile (34:35));
DD    = str2num (mms_ql__EDI__BdvE__dataFile (36:37));

%{
				 0         0         0         0         0
mms2_edi_slow_ql_efield_20150509_v0.1.4.cdf
	'Epoch'                          [1x2 double] [ 501] 'tt2000' 'T/'  'Full' 'None' [0] [-9223372036854775808]
	'Epoch_delta_plus'               [1x2 double] [   1] 'tt2000' 'F/'  'Full' 'None' [0] [-9223372036854775808]
	'epoch_gd12_beam'                [1x2 double] [3488] 'tt2000' 'T/'  'Full' 'None' [0] [-9223372036854775808]
	'epoch_gd21_beam'                [1x2 double] [3142] 'tt2000' 'T/'  'Full' 'None' [0] [-9223372036854775808]
	'mms2_edi_E_dmpa'                [1x2 double] [ 501] 'single' 'T/T' 'Full' 'None' [0] [      -1.0000000e+30]
	'mms2_edi_v_ExB_dmpa'            [1x2 double] [ 501] 'single' 'T/T' 'Full' 'None' [0] [      -1.0000000e+30]
	'mms2_edi_d_dmpa'                [1x2 double] [ 501] 'single' 'T/T' 'Full' 'None' [0] [      -1.0000000e+30]
	'mms2_edi_B_dmpa'                [1x2 double] [ 501] 'single' 'T/T' 'Full' 'None' [0] [      -1.0000000e+30]
	'mms2_edi_quality'               [1x2 double] [ 501] 'int8'   'T/'  'Full' 'None' [0] [                -127]
	'mms2_edi_pos_virtual_gun1_dmpa' [1x2 double] [3488] 'single' 'T/T' 'Full' 'None' [0] [      -1.0000000e+30]
	'mms2_edi_pos_virtual_gun2_dmpa' [1x2 double] [3142] 'single' 'T/T' 'Full' 'None' [0] [      -1.0000000e+30]
	'mms2_edi_fv_gd12_dmpa'          [1x2 double] [3488] 'single' 'T/T' 'Full' 'None' [0] [      -1.0000000e+30]
	'mms2_edi_fv_gd21_dmpa'          [1x2 double] [3142] 'single' 'T/T' 'Full' 'None' [0] [      -1.0000000e+30]
	'mms2_edi_recnum'                [1x2 double] [ 501] 'int32'  'T/'  'Full' 'None' [0] [         -2147483647]
	'mms2_edi_recnum_gd12'           [1x2 double] [3488] 'int32'  'T/'  'Full' 'None' [0] [         -2147483647]
	'mms2_edi_recnum_gd21'           [1x2 double] [3142] 'int32'  'T/'  'Full' 'None' [0] [         -2147483647]
	'mms2_edi_beam_quality_gd12'     [1x2 double] [3488] 'int8'   'T/'  'Full' 'None' [0] [                -127]
	'mms2_edi_beam_quality_gd21'     [1x2 double] [3142] 'int8'   'T/'  'Full' 'None' [0] [                -127]
	'mms2_edi_d_std_dmpa'            [1x2 double] [ 501] 'single' 'T/T' 'Full' 'None' [0] [      -1.0000000e+30]
	'mms2_edi_B_std_dmpa'            [1x2 double] [ 501] 'single' 'T/T' 'Full' 'None' [0] [      -1.0000000e+30]
	'E_Labl_Ptr'                     [1x2 double] [   3] 'char'   'T/'  'Full' 'None' [0] '  '
	'v_Labl_Ptr'                     [1x2 double] [   3] 'char'   'T/'  'Full' 'None' [0] '  '
	'd_Labl_Ptr'                     [1x2 double] [   3] 'char'   'T/'  'Full' 'None' [0] '  '
	'B_Labl_Ptr'                     [1x2 double] [   3] 'char'   'T/'  'Full' 'None' [0] '  '
	'vg_labl_vname'                  [1x2 double] [   3] 'char'   'T/'  'Full' 'None' [0] ' '

				 0         0         0         0         0
mms2_edp_comm_ql_dce2d_20150509120000_v0.1.0.cdf
  'mms2_edp_dce_epoch'      [1x2 double] [2764735] 'tt2000' 'T/'  'Full' 'None'   [    0] [-9223372036854775808]
  'LABL_1'                  [1x2 double] [      1] 'char'   'F/T' 'Full' 'GZIP.6' [    1] '     '
  'mms2_edp_dce_xyz_dsl'    [1x2 double] [2764735] 'single' 'T/T' 'Full' 'GZIP.6' [ 5462] [      -1.0000000e+30]
  'mms2_edp_dce_bitmask'    [1x2 double] [2764735] 'uint8'  'T/'  'Full' 'GZIP.6' [65536] [                 254]
  'mms2_edp_dce_quality'    [1x2 double] [2764735] 'int16'  'T/'  'Full' 'GZIP.6' [32768] [              -32767]
%}

% ~~~~~~~~~~~~~~~~~~~
% If CDF_FileInfo.FileSettings.Majority = 'Row' then transpose to col vector matrix.
% But MATLAB searches faster in cols, so if there is filtering to be done,
% better to do it along cols, and transpose only when necessary for math.
% ~~~~~~~~~~~~~~~~~~~

if ~isequal (mms_ql__EDI__BdvE__dataFile, 0) % then a valid [we hope] file selected
	disp ([ 'Reading MMS EDI_&_B data... ', mms_ql__EDI__BdvE__data ])
	mms_ql_EDI_dataFile_info = spdfcdfinfo (mms_ql__EDI__BdvE__data);
	edi_B_sdcs_varInfo = CDF_varInfo (mms_ql_EDI_dataFile_info, ['mms', obsID, '_edi_B_dmpa']);
	edi_B_sdcs_FillVal = edi_B_sdcs_varInfo.FillVal;

	% ~~~~~~~~~~~~~~~~~~~ B d v E
	edi_BdvE_t2k = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              'Epoch', ...
		'ConvertEpochToDatenum', false, ...
		'KeepEpochAsIs',         true);
	edi_B_sdcs = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_B_dmpa']);
	edi_E_sdcs = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_E_dmpa']);
	edi_E_sdcs = edi_E_sdcs';

	% The record number of edi_B_sdcs to which 5s sets of EDI beams map.
	edi_BdvE_recnum = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_recnum']);
	disp 'Date range of edi_BdvE_t2k'
	[ datestr(spdftt2000todatenum(edi_BdvE_t2k(1)),   'yyyy-mm-dd HH:MM:ss'), ' ',...
	  datestr(spdftt2000todatenum(edi_BdvE_t2k(end)), 'yyyy-mm-dd HH:MM:ss') ]

	% ~~~~~~~~~~~~~~~~~~~ GDU data
	edi_gd12_beam_t2k = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              'epoch_gd12_beam', ...
		'ConvertEpochToDatenum', false, ...
		'KeepEpochAsIs',         true);
	edi_gd12_virtual_sdcs = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_pos_virtual_gun1_dmpa']);
	edi_gd12_fv_sdcs = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_fv_gd12_dmpa']);
	edi_gd12_xref2_B = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_recnum_gd12']);
% 	edi_gd12_quality = spdfcdfread (mms_ql__EDI__BdvE__data, ...
% 		'CombineRecords',        true, ...
% 		'Variable',              ['mms', obsID, '_edi_beam_quality_gd12']);
	disp 'Date range of edi_gd12_beam_t2k'
	[ datestr(spdftt2000todatenum(edi_gd12_beam_t2k(1)),   'yyyy-mm-dd HH:MM:ss'), ' ',...
	  datestr(spdftt2000todatenum(edi_gd12_beam_t2k(end)), 'yyyy-mm-dd HH:MM:ss') ]

	% ~~~~~~~~~~~~~~~~~~~ GDU data
	edi_gd21_beam_t2k = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              'epoch_gd21_beam', ...
		'ConvertEpochToDatenum', false, ...
		'KeepEpochAsIs',         true);
	edi_gd21_virtual_sdcs = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_pos_virtual_gun2_dmpa']);
	edi_gd21_fv_sdcs = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_fv_gd21_dmpa']);
	edi_gd21_xref2_B = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_recnum_gd21']);
% 	edi_gd21_quality = spdfcdfread (mms_ql__EDI__BdvE__data, ...
% 		'CombineRecords',        true, ...
% 		'Variable',              ['mms', obsID, '_edi_beam_quality_gd21']);
	disp 'Date range of edi_gd21_beam_t2k'
	[ datestr(spdftt2000todatenum(edi_gd21_beam_t2k(1)),   'yyyy-mm-dd HH:MM:ss'), ' ',...
	  datestr(spdftt2000todatenum(edi_gd21_beam_t2k(end)), 'yyyy-mm-dd HH:MM:ss') ]

	% If CDF_FileInfo.FileSettings.Majority = 'Row' then transpose to col vector matrix.
	% But MATLAB searches faster in cols, so if there is filtering to be done,
	% better to do it along cols, and transpose only when necessary for math.

	% ~~~~~~~~~~~~~~~~~~~
	iB_eq_NaN = find (isnan (edi_B_sdcs (:,1)));
	if ~isempty (iB_eq_NaN)
		iB_eq_NaN % debug CDF file for NaNs.
	end

	iBeqFillVal = find (edi_B_sdcs == edi_B_sdcs_FillVal);
	iBeqBad = union (iB_eq_NaN, iBeqFillVal);
	% It's OK to delete these 'bad' BdvE records, because edi_BdvE_recnum keeps track of the remainder
	edi_BdvE_t2k (iBeqBad   )    = [];
	edi_B_sdcs      (iBeqBad, :) = [];
	edi_BdvE_recnum (iBeqBad, :) = [];
	% At this point, all BdvE, tt2000, edi_BdvE_recnum records are in sync

	% Now is the time to change from nx3 data to 3xn.
	edi_gd_beam_t2k     = [ edi_gd12_beam_t2k'        edi_gd21_beam_t2k' ];
	edi_gd_virtual_sdcs = [ edi_gd12_virtual_sdcs'    edi_gd21_virtual_sdcs' ];
	edi_gd_fv_sdcs      = [ double(edi_gd12_fv_sdcs') double(edi_gd21_fv_sdcs') ];
	% edi_xref2_BdvE[] values match values in edi_BdvE_recnum[].
	% To match BdvE records with EDI records (1:many), either
	% 1) find (edi_BdvE_recnum == edi_xref2_BdvE(i)) ~> should return only 1 record. or
	% 2) find (edi_xref2_BdvE == edi_BdvE_recnum(i)) ~> can return 0 or more
	edi_xref2_BdvE      = [ edi_gd12_xref2_B'         edi_gd21_xref2_B' ];

	% Better keep track of which GDU fired the beam. Corresponds to concatenation order above.
	edi_gd_ID           = [ zeros(1, size(edi_gd12_fv_sdcs, 1), 'uint8')+1, zeros(1, size(edi_gd21_fv_sdcs, 1), 'uint8')+2 ];

	% ... and tranpose B records
	BdvE_dn    = spdftt2000todatenum (edi_BdvE_t2k);
	edi_BdvE_t2k = edi_BdvE_t2k';
	edi_B_sdcs      = edi_B_sdcs';

	[ ~, iSorted_beam_t2k ] = sort (edi_gd_beam_t2k, 2);
	% 	[ edi_gd_beam_t2k(iSorted_beam_t2k(1:20))', ...
	% 	  edi_gd_ID(iSorted_beam_t2k(1:20))'        ]

	% ~~~~~~~~~~~~~~~~~~~
	% mms2_edp_comm_ql_dce2d_20150509120000_v0.1.0.cdf'
	% ~~~~~~~~~~~~~~~~~~~
	if UseFileOpenGUI
		[mms_ql__EDP_dataFile, mms_ql_dataPath] = uigetfile ('mms*.cdf', 'Select an MMS ql EDP CDF file');
		if isequal (mms_ql__EDP_dataFile,  0) % then no valid file selected
			msgbox ('No valid MMS ql EDP data file selected.');
		else
			mms_ql__EDP_data = [mms_ql_dataPath, mms_ql__EDP_dataFile];
		end
	end

	if ~isequal (mms_ql__EDP_dataFile, 0) % then a valid [we hope] file selected
		disp ([ 'Reading EDP data... ', mms_ql__EDP_data ])
		mms_ql_EDP_dataFile_info = spdfcdfinfo (mms_ql__EDP_data);
		edp_dce_varInfo = CDF_varInfo (mms_ql_EDP_dataFile_info, ['mms', obsID, '_edp_dce_xyz_dsl']);
		edp_dce_fillVal = edp_dce_varInfo.FillVal;
% keyboard
		% ~~~~~~~~~~~~~~~~~~~ DC E-field
		edp_t2k = spdfcdfread (mms_ql__EDP_data, ...
			'CombineRecords',        true, ...
			'Variable',              'mms2_edp_dce_epoch', ...
			'ConvertEpochToDatenum', false, ...
			'KeepEpochAsIs',         true);
		% Electric field in DSL coordinates (DSL ~= DMPA ~= DBCS)
		edp_dce_xyz_sdcs = spdfcdfread (mms_ql__EDP_data, ...
			'CombineRecords',        true, ...
			'Variable',              ['mms', obsID, '_edp_dce_xyz_dsl']);

		idceFillVal = find (edp_dce_xyz_sdcs (:, 1) < edp_dce_fillVal);
		edp_dce_xyz_sdcs (idceFillVal, :) = [];
		edp_t2k      (idceFillVal)    = [];
		edp_dn = spdftt2000todatenum (edp_t2k);
% 		disp 'Date range of edp_t2k'
% 		[ datestr(spdftt2000todatenum(edp_t2k(1)),   'yyyy-mm-dd HH:MM:ss'), ' ',...
% 		  datestr(spdftt2000todatenum(edp_t2k(end)), 'yyyy-mm-dd HH:MM:ss') ]

		% keep EDP data that is in the range of EDI data
		iEDP_lt_EDI = find (edp_dn < BdvE_dn (1));
		iEDP_gt_EDI = find (edp_dn > BdvE_dn (end));
		iEDP_EDI_noMatch = union (iEDP_lt_EDI, iEDP_gt_EDI);
		edp_dn (iEDP_EDI_noMatch) = [];
		edp_dce_xyz_sdcs (iEDP_EDI_noMatch, :) = [];
		disp 'Date range of edp_t2k'
		[ datestr(spdftt2000todatenum(edp_t2k(1)),   'yyyy-mm-dd HH:MM:ss'), ' ',...
		  datestr(spdftt2000todatenum(edp_t2k(end)), 'yyyy-mm-dd HH:MM:ss') ]

		ValidDataLoaded = true; % assume that user knows data is good
	end

end % ~isequal (mms_ql__EDI__BdvE__dataFile, 0)
