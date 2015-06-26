% edi_drift_step_read_ql__EDI__B__data

UseFileOpenGUI = true;
if UseFileOpenGUI
	[mms_ql__EDI__B__dataFile, mms_ql_dataPath] = uigetfile ('mms*.cdf', 'Select an MMS ql EDI_&_B CDF file');
	if isequal (mms_ql__EDI__B__dataFile,  0) % then no valid file selected
		msgbox ('No valid CSA EFW L2 data file selected.');
	else
		mms_ql__EDI__B__data = [mms_ql_dataPath, mms_ql__EDI__B__dataFile];
	end
end

%{
	         0         0         0         0         0
	mms2_edi_slow_ql_efield-5min_20150509160800_v0.0.0.cdf
	2015-05-09 14:10:15 to 2015-05-09 16:31:40
	'Epoch'                      [1x2 double] [ 60] 'tt2000' 'T/' 'Full' 'None' [0] [-9223372036854775808]
	'epoch_gd12_beam'            [1x2 double] [682] 'tt2000' 'T/' 'Full' 'None' [0] [-9223372036854775808]
	'epoch_gd21_beam'            [1x2 double] [668] 'tt2000' 'T/' 'Full' 'None' [0] [-9223372036854775808]
	'mms*_edi_E_dmpa'            [1x2 double] [ 60] 'single' 'T/T 'Full' 'None' [0] [      -1.0000000e+30]
	'mms*_edi_v_ExB_dmpa'        [1x2 double] [ 60] 'single' 'T/T 'Full' 'None' [0] [      -1.0000000e+30]
	'mms*_edi_d_dmpa'            [1x2 double] [ 60] 'single' 'T/T 'Full' 'None' [0] [      -1.0000000e+30]
	'mms*_edi_B_dmpa'            [1x2 double] [ 60] 'single' 'T/T 'Full' 'None' [0] [      -1.0000000e+30]
	'mms*_edi_fv_gd12_dmpa'      [1x2 double] [682] 'single' 'T/T 'Full' 'None' [0] [      -1.0000000e+30]
	'mms*_edi_fv_gd21_dmpa'      [1x2 double] [668] 'single' 'T/T 'Full' 'None' [0] [      -1.0000000e+30]
	'mms*_edi_inds_gd12'         [1x2 double] [682] 'int32'  'T/' 'Full' 'None' [0] [         -2147483647]
	'mms*_edi_inds_gd21'         [1x2 double] [668] 'int32'  'T/' 'Full' 'None' [0] [         -2147483647]
	'mms*_edi_beam_quality_gd12' [1x2 double] [682] 'int32'  'T/' 'Full' 'None' [0] [         -2147483647]
	'mms*_edi_beam_quality_gd21' [1x2 double] [668] 'int32'  'T/' 'Full' 'None' [0] [         -2147483647]
	'E_Labl_Ptr'                 [1x2 double] [  3] 'char'   'T/' 'Full' 'None' [0] '  '
	'v_Labl_Ptr'                 [1x2 double] [  3] 'char'   'T/' 'Full' 'None' [0] '  '
	'd_Labl_Ptr'                 [1x2 double] [  3] 'char'   'T/' 'Full' 'None' [0] '  '
	'B_Labl_Ptr'                 [1x2 double] [  3] 'char'   'T/' 'Full' 'None' [0] '  '

  'Epoch'                             [1x2 double]    [ 743]    'tt2000'    'T/'     'Full'    'None'    [0]    [-9223372036854775808]
  'epoch_gd12_beam'                   [1x2 double]    [3488]    'tt2000'    'T/'     'Full'    'None'    [0]    [-9223372036854775808]
  'epoch_gd21_beam'                   [1x2 double]    [3142]    'tt2000'    'T/'     'Full'    'None'    [0]    [-9223372036854775808]
  'mms2_edi_E_dmpa'                   [1x2 double]    [ 743]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'mms2_edi_v_ExB_dmpa'               [1x2 double]    [ 743]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'mms2_edi_d_dmpa'                   [1x2 double]    [ 743]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'mms2_edi_B_dmpa'                   [1x2 double]    [ 743]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'mms2_edi_pos_virtual_gun1_dmpa'    [1x2 double]    [3488]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'mms2_edi_pos_virtual_gun2_dmpa'    [1x2 double]    [3142]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'mms2_edi_fv_gd12_dmpa'             [1x2 double]    [3488]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'mms2_edi_fv_gd21_dmpa'             [1x2 double]    [3142]    'single'    'T/T'    'Full'    'None'    [0]    [      -1.0000000e+30]
  'mms2_edi_inds_gd12'                [1x2 double]    [3488]    'int32'     'T/'     'Full'    'None'    [0]    [         -2147483647]
  'mms2_edi_inds_gd21'                [1x2 double]    [3142]    'int32'     'T/'     'Full'    'None'    [0]    [         -2147483647]
  'mms2_edi_beam_quality_gd12'        [1x2 double]    [3488]    'int8'      'T/'     'Full'    'None'    [0]    [                -127]
  'mms2_edi_beam_quality_gd21'        [1x2 double]    [3142]    'int8'      'T/'     'Full'    'None'    [0]    [                -127]
  'E_Labl_Ptr'                        [1x2 double]    [   3]    'char'      'T/'     'Full'    'None'    [0]    '  '                  
  'v_Labl_Ptr'                        [1x2 double]    [   3]    'char'      'T/'     'Full'    'None'    [0]    '  '                  
  'd_Labl_Ptr'                        [1x2 double]    [   3]    'char'      'T/'     'Full'    'None'    [0]    '  '                  
  'B_Labl_Ptr'                        [1x2 double]    [   3]    'char'      'T/'     'Full'    'None'    [0]    '  '                  
  'vg_labl_vname'                     [1x2 double]    [   3]    'char'      'T/'     'Full'    'None'    [0]    ' '                   
%}

obsID = mms_ql__EDI__B__dataFile (4:4);
YYYY  = str2num (mms_ql__EDI__B__dataFile (30:33));
MM    = str2num (mms_ql__EDI__B__dataFile (34:35));
DD    = str2num (mms_ql__EDI__B__dataFile (36:37));

if ~isequal (mms_ql__EDI__B__dataFile,  0) % then a valid [we hope] file selected
	disp ([ 'Reading MMS EDI_&_B data... ', mms_ql__EDI__B__data ])
	mms_ql_dataFile_info = spdfcdfinfo (mms_ql__EDI__B__data);

	edi_BdvE_tt2000 = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              'Epoch', ...
		'ConvertEpochToDatenum', false, ...
		'KeepEpochAsIs',         true);
	edi_B_dmpa = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_B_dmpa']);
	[ datestr(spdftt2000todatenum(edi_BdvE_tt2000(1)),   'yyyy-mm-dd HH:MM:ss'), ' ',...
	  datestr(spdftt2000todatenum(edi_BdvE_tt2000(end)), 'yyyy-mm-dd HH:MM:ss') ]

	edi_gd12_beam_tt2000 = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              'epoch_gd12_beam', ...
		'ConvertEpochToDatenum', false, ...
		'KeepEpochAsIs',         true);
	edi_gd12_virtual_dmpa = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_pos_virtual_gun1_dmpa']);
	edi_gd12_fv_dmpa = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_fv_gd12_dmpa']);
	edi_gd12_B_index = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_inds_gd12']);
	edi_gd12_quality = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_beam_quality_gd12']);

	edi_gd21_beam_tt2000 = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              'epoch_gd21_beam', ...
		'ConvertEpochToDatenum', false, ...
		'KeepEpochAsIs',         true);
	edi_gd21_virtual_dmpa = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_pos_virtual_gun2_dmpa']);
	edi_gd21_fv_dmpa = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_fv_gd21_dmpa']);
	edi_gd21_B_index = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_inds_gd21']);
	edi_gd21_quality = spdfcdfread (mms_ql__EDI__B__data, ...
		'CombineRecords',        true, ...
		'Variable',              ['mms', obsID, '_edi_beam_quality_gd21']);

	% If CDF_FileInfo.FileSettings.Majority = 'Row' then transpose to col vector matrix.
	% But MATLAB searches faster in cols, so if there is filtering to be done,
	% better to do it along cols, and transpose only at when necessary for math.

	% We want only quality=3 beams - rumor is that the file will have only Q=3 data
	% The data arrives in row vectors, so index via (index, :)
% 	iq12lt3 = find (edi_gd12_quality < 3);
% 	iq21lt3 = find (edi_gd21_quality < 3);
	% Now do the same for NaNs?

% 	if ~isempty (iq12lt3)
% 		size (edi_gd12_beam_tt2000)
% 		edi_gd12_beam_tt2000 (iq12lt3, :) = [];
% 		edi_gd12_fv_dmpa     (iq12lt3, :) = [];
% 		edi_gd12_B_index     (iq12lt3, :) = [];
% 		edi_gd12_quality     (iq12lt3, :) = [];
% 		size (edi_gd12_beam_tt2000)
% 	end

% 	if ~isempty (iq21lt3)
% 		size (edi_gd21_beam_tt2000)
% 		edi_gd21_beam_tt2000 (iq21lt3, :) = [];
% 		edi_gd21_fv_dmpa     (iq21lt3, :) = [];
% 		edi_gd21_B_index     (iq21lt3, :) = [];
% 		edi_gd21_quality     (iq21lt3, :) = [];
% 		size (edi_gd21_beam_tt2000)
% 	end

	iBeqNaN = find (isnan (edi_B_dmpa (:,1)));
	if ~isempty (iBeqNaN)
		iBeqNaN
	end
	edi_B_dmpa_FillVal = -1.0000000e+30;
	iBeqFillVal = find (edi_B_dmpa == edi_B_dmpa_FillVal);
	iBeqBad = union (iBeqNaN, iBeqFillVal);
	edi_BdvE_tt2000 (iBeqBad   ) = [];
	edi_B_dmpa     (iBeqBad, :) = [];

	% There may be NaNs in the B data
	% If any 5 s period has no B data, we can't use those BEAMs as CENTER beams,
	% but we could use them as beams on either side, so the logic gets more complicated...
	% If we remove the bad B data, we must remove the EDI data that corresponds.
	% If we keep the bad B records, then we must check B valid for each EDI record set.
	% The problem with deleting the bad B data records is that we mess up the 
	% EDI-record:B-record index is corrupted, because there is NO B-record index ---
	% the EDI index cross references are based on B-record position. Therefore,
	% if we delete a B-record, all the following records get bumped down by one,
	% and this kills the EDI:B record cross-reference.

	clear edi_gd12_quality
	clear edi_gd21_quality

	% Now the problem is that we want to scroll thru the firing vectors (beams)
	% in linear time, choosing the center beam, and finding 4-8 beams on either side
	% FROM EITHER GDU.
	% That means concat time, fv, and B_index. Q is no longer necessary.
	% Then create index of time, sorted in ascending order. Scroll thru the data
	% during analysis, using the sorted time index.

	% Now is the time to change from nx3 data to 3xn.
	edi_gd_beam_tt2000  = [ edi_gd12_beam_tt2000'     edi_gd21_beam_tt2000' ];
	edi_gd_virtual_dmpa = [ edi_gd12_virtual_dmpa'    edi_gd21_virtual_dmpa' ];
	edi_gd_fv_dmpa      = [ double(edi_gd12_fv_dmpa') double(edi_gd21_fv_dmpa') ];
	edi_gd_B_index      = [ edi_gd12_B_index'         edi_gd21_B_index' ];
	% Better keep track of which GDU fired the beam. Corresponds to concatenation order above.
	edi_gd_ID           = [ zeros(1,length(edi_gd12_fv_dmpa),'uint8')+1, zeros(1,length(edi_gd21_fv_dmpa),'uint8')+2 ];
	% ... and tranpose B records
	edi_BdvE_tt2000 = edi_BdvE_tt2000';
	edi_B_dmpa      = edi_B_dmpa';

	[ ~, iSorted_beam_tt2000 ] = sort (edi_gd_beam_tt2000, 2);
	% 	[ edi_gd_beam_tt2000(iSorted_beam_tt2000(1:20))', ...
	% 	  edi_gd_ID(iSorted_beam_tt2000(1:20))'        ]
	ValidDataLoaded = true; % assume that user knows data is good

end
