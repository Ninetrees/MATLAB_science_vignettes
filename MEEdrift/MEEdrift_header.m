%{
	Name: MEEdrift
	Description: Correct EFW E-field measurements based on EDI measurements
	Revision History (Version; Author; Date; Reason):
	v0100; 2015-06-16;
		V 01.00
		Uses special EDI, B, EDP extracts from commissioning data for initial dev
		MMS2: May 9, 2015  16:08 - 16:13
		MMS4: May 6, 2015  15:30 - 15:35
		Prepare MEEdrift for git, GitHub version control system

	Important !!! notes about variable naming conventions
	iVarName is an index to VarName
	iiVarName is an index to an index to VarName; e.g., VarName (iVarName (iiVarName))
	nVarName is the number of elements in VarName
	_dn  indicates a MATLAB datenum variable
	_t2k indicates a CDF TT2000 variable
	u or _u is a unit vector; _u used for clarity when needed

	MATLAB datenum does not include leap seconds as of 2015-07-28, so...
	the chioces are tt2000 or tt2000~>ssm (seconds since midnight).
	ssm has the advantage of already being in seconds, where so many calcs are done,
	but doesn't cross midnight without handling at each calc.
	tt2000 handles leap seconds and midnight, but requires handling in some cases;
	e.g., if abs ((edp_t2k (iEDP_t2kHI) - edi_gd_beam_t2k (iBeam)) < 1.0e8),
	rather than if abs ((edp_ssm (iEDP_t2kHI) - edi_gd_beam_ssm (iBeam)) < 0.1s).
	In this cose, note the uniform conversion from tt2k in ns to some other time
	simply by using the relationship tt2k * 1e9 = seconds.

	A note about MATLAB & matrices in MEEdrift:
	!!! Matrices are row vectors for efficiency 'find'ing, plotting !!!
	!!! but all math is with col vectors for linear algebra consistency !!!
	!!! This means that when vectors are filled from row vector matrices, they are transposed. !!!
	In MATLAB, a vector is a one-dimensional array and a matrix is a two-dimensional array.
	MDAa are created by extending 2D arrays; e.g., [:, :, n, ...].
		For example: (3, 2, 5, nEDP) =
			a 3D vector, for each of 2 intersecting beams, for up to 5 matches, for each EFW event
	MEEdrift 1D row or col vectors are addressed simply by V (i).
	MEEdrift 2D & 3D vectors are considered column vectors; e.g., V = [ v1 v2 ...], where each vector is v1 = [ v11; v12; ... ],
	unless noted.
	MEEdrift rotation matrices are 3x3 matrices of 3D vectors; arrays of them are 3x3xn.
	spdfcdfread places data into row vectors, so matrices returned as n x 2 (for example), but we typically use 2 x n.

	For efficiency, matrices are pre-allocated using syntax similar to: M = zeros (1, n_records, 'double'),
	UNLESS routines such as textread or spdfcdfread allocate space during the read process.

	For clarification, variables are cleared when no longer needed.
	Generally, variables beginning with 'i' are indices, 'n' are counts.

	Assume
	~ 1° degree knowledge of beam firing direction
	~ Time of flight of the electrons with a resolution better than 1 us

	Viewed from +GSEz, the MMS spacecraft rotate CCW, and spacecraft coords are GSE-like;
	i.e., DBCS ~= DMPA ~= DSL, and all 3 are close to GSE.

	!!! CDF record numbers are zer0-based !!! MATLAB records numbers are 1-based. Adjust accordingly.

	mms2_edi_slow_ql_efield_20150509_v0.1.4.cdf implements a 1:many record cross-reference twixt
	mms2_edi_recnum:mms2_edi_recnum_gd12|mms2_edi_recnum_gd21. These xrefs are used to link EDI records to
	B(dvE) records.
%}
