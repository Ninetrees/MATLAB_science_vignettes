% MEEdrift_read_ql__EDI__B__EDP_data

myLibCDFConstants

% ~~~~~~~~~~~~~~~~~~~
% mms2_edi_slow_ql_efield_20150509_v0.1.4.cdf'
% ~~~~~~~~~~~~~~~~~~~
if UseFileOpenGUI
	[mms_ql__EDI__BdvE__dataFile, mms_ql_dataPath] = uigetfile ('mms*.cdf', 'Select an MMS ql EDI_&_B CDF file');
	if isequal (mms_ql__EDI__BdvE__dataFile,  0) % then no valid file selected
		msgbox ('No valid MMS ql EDI data file selected.');
	else
		mms_ql__EDI__BdvE__data = [mms_ql_dataPath cFileSep mms_ql__EDI__BdvE__dataFile];
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
mms2_edp_comm_ql_dce2d_20150509120000_v0.1.0.cdf'
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
	mms_ql_EDI_fileInfo = spdfcdfinfo (mms_ql__EDI__BdvE__data);
	edi_B_dmpa_varInfo = CDF_varInfo (mms_ql_EDI_fileInfo, ['mms', obsID, '_edi_B_dmpa']);
	edi_B_dmpa_fillVal = edi_B_dmpa_varInfo.fillVal;
	edi_d_dmpa_varInfo = CDF_varInfo (mms_ql_EDI_fileInfo, ['mms', obsID, '_edi_d_dmpa']);
	edi_d_dmpa_fillVal = edi_d_dmpa_varInfo.fillVal;

	% ~~~~~~~~~~~~~~~~~~~ B d v E
	edi_BdvE_t2k = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              'Epoch', ...
		'ConvertEpochToDatenum', false, ...
		'KeepEpochAsIs',         true);
	edi_B_dmpa = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_B_dmpa']);
	edi_E_dmpa = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_E_dmpa']);
	edi_d_dmpa = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_d_dmpa']);

	% The record number of edi_B_dmpa to which 5s sets of EDI beams map.
	BdvE_xref2_edi = spdfcdfread (mms_ql__EDI__BdvE__data, ...
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
	edi_gd12_virtual_dmpa = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_pos_virtual_gun1_dmpa']);
	edi_gd12_fv_dmpa = spdfcdfread (mms_ql__EDI__BdvE__data, ...
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
	edi_gd21_virtual_dmpa = spdfcdfread (mms_ql__EDI__BdvE__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_pos_virtual_gun2_dmpa']);
	edi_gd21_fv_dmpa = spdfcdfread (mms_ql__EDI__BdvE__data, ...
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
	iB_eq_NaN = find (isnan (edi_B_dmpa (:,1)));
	if ~isempty (iB_eq_NaN)
		iB_eq_NaN % debug CDF file for NaNs.
	end

	iBeqFillVal = find (edi_B_dmpa == edi_B_dmpa_fillVal);
	iBeqBad = union (iB_eq_NaN, iBeqFillVal);
	% It's OK to delete these 'bad' BdvE records, because edi_BdvE_recnums keeps track of the remainder
	edi_BdvE_t2k     (iBeqBad   ) = [];
	edi_B_dmpa       (iBeqBad, :) = [];
	edi_BdvE_recnums (iBeqBad, :) = [];
	% At this point, all BdvE, tt2000, edi_BdvE_recnums records are in sync

	% There may be NaNs or FillVals in the B data
	% If any 5 s period has no B data, we can't use those BEAMs as CENTER beams,
	% but we could use them as beams on either side.
	% If we remove the bad B data, we must remove the EDI data that corresponds.
	% If we keep the NaN B records, then we must check B valid for each EDI record set.

	for i = 1: length (iBeqBad)
		i12BeqBad_match_EDI = find (edi_gd12_xref2_B == iBeqBad (i));
		if ~isempty (i12BeqBad_match_EDI)
			edi_gd12_beam_t2k     (i12BeqBad_match_EDI, :) = [];
			edi_gd12_fv_dmpa      (i12BeqBad_match_EDI, :) = [];
			edi_gd12_virtual_dmpa (i12BeqBad_match_EDI, :) = [];
			edi_gd12_xref2_B      (i12BeqBad_match_EDI, :) = [];  %%% See note above
		end

		i21BeqBad_match_EDI = find (edi_gd21_xref2_B == iBeqBad (i));
		if ~isempty (i21BeqBad_match_EDI)
			edi_gd21_beam_t2k     (i21BeqBad_match_EDI, :) = [];
			edi_gd21_fv_dmpa      (i21BeqBad_match_EDI, :) = [];
			edi_gd21_virtual_dmpa (i21BeqBad_match_EDI, :) = [];
			edi_gd21_xref2_B      (i21BeqBad_match_EDI, :) = [];  %%% See note above
		end
	end

	% We want to scroll thru the firing vectors (beams) in linear time,
	% choosing the center beam, and finding 4-8 beams on either side
	% FROM EITHER GDU.
	% That means concat time, fv, and B_index. Q is no longer necessary.
	% Then create index of time, sorted in ascending order. Scroll thru the data
	% during analysis, using the sorted time index.

% 	clear edi_gd12_quality % we use it only to filter initially
% 	clear edi_gd21_quality

	% Now is the time to change from nx3 data to 3xn, where applicable
	edi_gd_beam_t2k = [ edi_gd12_beam_t2k; edi_gd21_beam_t2k ];
	edi_gd_beam_dn  = spdftt2000todatenum (edi_gd_beam_t2k);

	edi_gd_fv_dmpa      = [ double(edi_gd12_fv_dmpa)' double(edi_gd21_fv_dmpa)' ];
	edi_gd_virtual_dmpa = [ edi_gd12_virtual_dmpa'    edi_gd21_virtual_dmpa' ];
	% edi_xref2_BdvE[] values match values in edi_BdvE_recnums[].
	% To match BdvE records with EDI records (1:many), either
	% 1) find (edi_BdvE_recnums == edi_xref2_BdvE(i)) ~> should return only 1 record. or
	% 2) find (edi_xref2_BdvE == edi_BdvE_recnums(i)) ~> can return 0 or more
	edi_xref2_BdvE      = [ edi_gd12_xref2_B;         edi_gd21_xref2_B ];

	% Better keep track of which GDU fired the beam. Corresponds to concatenation order above.
	edi_gd_ID           = [ ...
		zeros(1, size(edi_gd12_fv_dmpa, 1), 'uint8')+1, ...
		zeros(1, size(edi_gd21_fv_dmpa, 1), 'uint8')+2 ];

	% ... and tranpose B records
	BdvE_dn    = spdftt2000todatenum (edi_BdvE_t2k);
	edi_B_dmpa = edi_B_dmpa';
	edi_d_dmpa = edi_d_dmpa';
	edi_E_dmpa = edi_E_dmpa';

	[ ~, iSorted_beam_t2k ] = sort (edi_gd_beam_t2k, 1);
	ValidDataLoaded = true; % assume that user knows data is good

	% ~~~~~~~~~~~~~~~~~~~
	% mms2_edp_comm_ql_dce2d_20150509120000_v0.1.0.cdf'
	% ~~~~~~~~~~~~~~~~~~~
	if UseFileOpenGUI
		[mms_ql__EDP__dataFile, mms_ql_dataPath] = uigetfile ('mms*.cdf', 'Select an MMS ql EDP CDF file');
		if isequal (mms_ql__EDP__dataFile,  0) % then no valid file selected
			msgbox ('No valid MMS ql EDP data file selected.');
		else
			mms_ql__EDP_data = [mms_ql_dataPath cFileSep mms_ql__EDP__dataFile];
		end
	end

	if ~isequal (mms_ql__EDP__dataFile, 0) % then a valid [we hope] file selected
		disp ([ 'Reading EDP data... ', mms_ql__EDP_data ])
		mms_ql_EDP_fileInfo = spdfcdfinfo (mms_ql__EDP_data);
		dce_xyz_dsl_varInfo = CDF_varInfo (mms_ql_EDP_fileInfo, ['mms', obsID, '_edp_dce_xyz_dsl']);
		dce_xyz_dsl_fillVal = dce_xyz_dsl_varInfo.fillVal;
% keyboard
		% ~~~~~~~~~~~~~~~~~~~ DC E-field
		edp_t2k = spdfcdfread (mms_ql__EDP_data, ...
			'CombineRecords',        true, ...
			'Variable',              'mms2_edp_dce_epoch', ...
			'ConvertEpochToDatenum', false, ...
			'KeepEpochAsIs',         true);
		% Electric field in DSL coordinates (DSL ~= DMPA ~= DBCS)
		edp_E3D_dsl = spdfcdfread (mms_ql__EDP_data, ...
			'CombineRecords',        true, ...
			'Variable',              ['mms', obsID, '_edp_dce_xyz_dsl']);
		edp_E3D_Q = spdfcdfread (mms_ql__EDP_data, ...
			'CombineRecords',        true, ...
			'Variable',              ['mms', obsID, '_edp_dce_quality']);

		idceFillVal = find (edp_E3D_dsl (:, 1) == dce_xyz_dsl_fillVal);
		edp_E3D_dsl (idceFillVal, :) = [];
		edp_t2k     (idceFillVal)    = [];
		edp_E3D_Q   (idceFillVal)    = [];

		iE3D_Q = find (edp_E3D_Q < 3); % 2015-07-30 they are all 0

		edp_dn = spdftt2000todatenum (edp_t2k);
% 		disp 'Date range of edp_t2k'
% 		[ datestr(spdftt2000todatenum(edp_t2k(1)),   'yyyy-mm-dd HH:MM:ss'), ' ',...
% 		  datestr(spdftt2000todatenum(edp_t2k(end)), 'yyyy-mm-dd HH:MM:ss') ]

		% keep EDP data that is in the range of EDI data
		iEDP_lt_EDI = find (edp_dn < BdvE_dn (1));
		iEDP_gt_EDI = find (edp_dn > BdvE_dn (end));
		iEDP_EDI_noMatch = union (iEDP_lt_EDI, iEDP_gt_EDI);

		edp_E3D_dsl (iEDP_EDI_noMatch, :) = [];
		edp_t2k     (iEDP_EDI_noMatch) = [];
		edp_E3D_Q   (iEDP_EDI_noMatch) = [];
		edp_dn      (iEDP_EDI_noMatch) = [];

		edp_E3D_dsl = edp_E3D_dsl';
% 		edp_E3D_Q   = edp_E3D_Q;
% 		edp_dn      = edp_dn;

		disp 'Date range of edp_t2k'
		[ datestr(spdftt2000todatenum(edp_t2k(1)),   'yyyy-mm-dd HH:MM:ss'), ' ',...
		  datestr(spdftt2000todatenum(edp_t2k(end)), 'yyyy-mm-dd HH:MM:ss') ]

		nBeams = size (edi_gd_fv_dmpa, 2);
		nEDP   = size (edp_E3D_dsl, 2);

		edp_zoomStart = 1;
		edp_zoomEnd   = 200;
		MEEdrift_plot_EDP_data

	end % if ~isequal (mms_ql__EDP__dataFile, 0)

end % ~isequal (mms_ql__EDI__BdvE__dataFile, 0)

	% We want only quality=3 beams; 2015-07-24: the data is only Q=3, so unnecessary
	% The data arrives in row vectors, so index via (index, :)
% 	iq12lt3 = find (edi_gd12_quality < 3);
% 	iq21lt3 = find (edi_gd21_quality < 3);
%
% 	if ~isempty (iq12lt3)
% 		size (edi_gd12_beam_t2k)
% 		edi_gd12_beam_t2k (iq12lt3, :) = [];
% 		edi_gd12_fv_dmpa     (iq12lt3, :) = [];
% 		edi_gd12_xref2_B     (iq12lt3, :) = [];
% 		edi_gd12_quality     (iq12lt3, :) = [];
% 		size (edi_gd12_beam_t2k)
% 	end
%
% 	if ~isempty (iq21lt3)
% 		size (edi_gd21_beam_t2k)
% 		edi_gd21_beam_t2k (iq21lt3, :) = [];
% 		edi_gd21_fv_dmpa     (iq21lt3, :) = [];
% 		edi_gd21_xref2_B     (iq21lt3, :) = [];
% 		edi_gd21_quality     (iq21lt3, :) = [];
% 		size (edi_gd21_beam_t2k)
% 	end
