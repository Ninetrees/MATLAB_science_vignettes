function varInfo = CDF_varInfo (fileInfo, varName)

	% CDF_fileInfo_.Variables
	% Example (2015-07-24): 'Epoch' [1x2 double] [ 501] 'tt2000' 'T/'  'Full' 'None' [0] [-9223372036854775808]

	% Not every fileInfo.VariableAttributes._attrib_ needs to exist, so relevant
	% warning ('CDF_varInfo: Variable "%s" _attrib_ not found in "%s".', varName, fileInfo.Filename)
	% are commented out.

	[nVariables, ~] = size (fileInfo.Variables);
	varFound = false;
	for iVariables = 1 : nVariables
		if strcmp (fileInfo.Variables {iVariables, 1}, varName)
			varInfo.Name        = fileInfo.Variables {iVariables, 1}; % use Caps?
			varInfo.Dims        = fileInfo.Variables {iVariables, 2};
			varInfo.nRecs       = fileInfo.Variables {iVariables, 3};
			varInfo.DataType    = fileInfo.Variables {iVariables, 4}; % DataType
			varInfo.DimsVary    = fileInfo.Variables {iVariables, 5}; % RecDimVariance
			varInfo.Composition = fileInfo.Variables {iVariables, 6};
			varInfo.Compression = fileInfo.Variables {iVariables, 7};
			varInfo.Blocking    = fileInfo.Variables {iVariables, 8};
			varInfo.PadVal      = fileInfo.Variables {iVariables, 9};

			[nCATDESC, ~] = size (fileInfo.VariableAttributes.CATDESC);
			CATDESC_found = false;
			for iCATDESC = 1 : nCATDESC
				if strcmp (fileInfo.VariableAttributes.CATDESC {iCATDESC, 1}, varName)
					varInfo.CatDesc  = fileInfo.VariableAttributes.CATDESC {iCATDESC, 2};
					CATDESC_found = true;
				end
			end
			if ~CATDESC_found
				varInfo.CatDesc = [];
				warning ('CDF_varInfo: Variable "%s" CATDESC not found in "%s".', varName, fileInfo.Filename)
			end

			[nDEPEND_0, ~] = size (fileInfo.VariableAttributes.DEPEND_0);
			DEPEND_0_found = false;
			for iDEPEND_0 = 1 : nDEPEND_0
				if strcmp (fileInfo.VariableAttributes.DEPEND_0 {iDEPEND_0, 1}, varName)
					varInfo.Depend_0  = fileInfo.VariableAttributes.DEPEND_0 {iDEPEND_0, 2};
					DEPEND_0_found = true;
				end
			end
			if ~DEPEND_0_found
				varInfo.Depend_0 = [];
				warning ('CDF_varInfo: Variable "%s" DEPEND_0 not found in "%s".', varName, fileInfo.Filename)
			end

% 			[nDEPEND_1, ~] = size (fileInfo.VariableAttributes.DEPEND_1);
% 			DEPEND_1_found = false;
% 			for iDEPEND_1 = 1 : nDEPEND_1
% 				if strcmp (fileInfo.VariableAttributes.DEPEND_1 {iDEPEND_1, 1}, varName)
% 					varInfo.Depend_1  = fileInfo.VariableAttributes.DEPEND_1 {iDEPEND_1, 2};
% 					DEPEND_1_found = true;
% 				end
% 			end
% 			if ~DEPEND_1_found
% 				varInfo.Depend_1 = [];
% 				warning ('CDF_varInfo: Variable "%s" DEPEND_1 not found in "%s".', varName, fileInfo.Filename)
% 			end

% 			[nDEPEND_2, ~] = size (fileInfo.VariableAttributes.DEPEND_2);
% 			DEPEND_2_found = false;
% 			for iDEPEND_2 = 1 : nDEPEND_2
% 				if strcmp (fileInfo.VariableAttributes.DEPEND_2 {iDEPEND_2, 1}, varName)
% 					varInfo.Depend_2  = fileInfo.VariableAttributes.DEPEND_2 {iDEPEND_2, 2};
% 					DEPEND_2_found = true;
% 				end
% 			end
% 			if ~DEPEND_2_found
% 				varInfo.Depend_2 = [];
% 				warning ('CDF_varInfo: Variable "%s" DEPEND_2 not found in "%s".', varName, fileInfo.Filename)
% 			end

% 			[nDEPEND_3, ~] = size (fileInfo.VariableAttributes.DEPEND_3);
% 			DEPEND_3_found = false;
% 			for iDEPEND_3 = 1 : nDEPEND_3
% 				if strcmp (fileInfo.VariableAttributes.DEPEND_3 {iDEPEND_3, 1}, varName)
% 					varInfo.Depend_3  = fileInfo.VariableAttributes.DEPEND_3 {iDEPEND_3, 2};
% 					DEPEND_3_found = true;
% 				end
% 			end
% 			if ~DEPEND_3_found
% 				varInfo.Depend_3 = [];
% 				warning ('CDF_varInfo: Variable "%s" DEPEND_3 not found in "%s".', varName, fileInfo.Filename)
% 			end

			[nFIELDNAM, ~] = size (fileInfo.VariableAttributes.FIELDNAM);
			FIELDNAM_found = false;
			for iFIELDNAM = 1 : nFIELDNAM
				if strcmp (fileInfo.VariableAttributes.FIELDNAM {iFIELDNAM, 1}, varName)
					varInfo.FieldName  = fileInfo.VariableAttributes.FIELDNAM {iFIELDNAM, 2};
					FIELDNAM_found = true;
				end
			end
			if ~FIELDNAM_found
				varInfo.FieldName = [];
				warning ('CDF_varInfo: Variable "%s" FIELDNAM not found in "%s".', varName, fileInfo.Filename)
			end

			[nFILLVAL, ~] = size (fileInfo.VariableAttributes.FILLVAL);
			FILLVAL_found = false;
			for iFILLVAL = 1 : nFILLVAL
				if strcmp (fileInfo.VariableAttributes.FILLVAL {iFILLVAL, 1}, varName)
					varInfo.FillVal  = fileInfo.VariableAttributes.FILLVAL {iFILLVAL, 2};
					FILLVAL_found = true;
				end
			end
			if ~FILLVAL_found
				varInfo.FillVal = [];
				warning ('CDF_varInfo: Variable "%s" FILLVAL not found in "%s".', varName, fileInfo.Filename)
			end

			[nLABLAXIS, ~] = size (fileInfo.VariableAttributes.LABLAXIS);
			LABLAXIS_found = false;
			for iLABLAXIS = 1 : nLABLAXIS
				if strcmp (fileInfo.VariableAttributes.LABLAXIS {iLABLAXIS, 1}, varName)
					varInfo.AxisLabel  = fileInfo.VariableAttributes.LABLAXIS {iLABLAXIS, 2};
					LABLAXIS_found = true;
				end
			end
			if ~LABLAXIS_found
				varInfo.AxisLabel = [];
% 				warning ('CDF_varInfo: Variable "%s" LABLAXIS not found in "%s".', varName, fileInfo.Filename)
			end

			[nSI_CONVERSION, ~] = size (fileInfo.VariableAttributes.SI_CONVERSION);
			SI_Conversion_found = false;
			for iSI_CONVERSION = 1 : nSI_CONVERSION
				if strcmp (fileInfo.VariableAttributes.SI_CONVERSION {iSI_CONVERSION, 1}, varName)
					varInfo.SI_conversion  = fileInfo.VariableAttributes.SI_CONVERSION {iSI_CONVERSION, 2};
					SI_Conversion_found = true;
				end
			end
			if ~SI_Conversion_found
				varInfo.SI_conversion = [];
				warning ('CDF_varInfo: Variable "%s" SI_CONVERSION not found in "%s".', varName, fileInfo.Filename)
			end

			[nUNITS, ~] = size (fileInfo.VariableAttributes.UNITS);
			Units_found = false;
			for iUNITS = 1 : nUNITS
				if strcmp (fileInfo.VariableAttributes.UNITS {iUNITS, 1}, varName)
					varInfo.Units  = fileInfo.VariableAttributes.UNITS {iUNITS, 2};
					Units_found = true;
				end
			end
			if ~Units_found
				varInfo.Units = [];
				warning ('CDF_varInfo: Variable "%s" UNITS not found in "%s".', varName, fileInfo.Filename)
			end

			[nVALIDMIN, ~] = size (fileInfo.VariableAttributes.VALIDMIN);
			Units_found = false;
			for iVALIDMIN = 1 : nVALIDMIN
				if strcmp (fileInfo.VariableAttributes.VALIDMIN {iVALIDMIN, 1}, varName)
					varInfo.ValidMin  = fileInfo.VariableAttributes.VALIDMIN {iVALIDMIN, 2};
					Units_found = true;
				end
			end
			if ~Units_found
				varInfo.ValidMin = [];
				warning ('CDF_varInfo: Variable "%s" VALIDMIN not found in "%s".', varName, fileInfo.Filename)
			end

			[nVALIDMAX, ~] = size (fileInfo.VariableAttributes.VALIDMAX);
			Units_found = false;
			for iVALIDMAX = 1 : nVALIDMAX
				if strcmp (fileInfo.VariableAttributes.VALIDMAX {iVALIDMAX, 1}, varName)
					varInfo.ValidMax  = fileInfo.VariableAttributes.VALIDMAX {iVALIDMAX, 2};
					Units_found = true;
				end
			end
			if ~Units_found
				varInfo.ValidMax = [];
				warning ('CDF_varInfo: Variable "%s" VALIDMAX not found in "%s".', varName, fileInfo.Filename)
			end

			varFound = true;
		end
	end
	if ~varFound
		varInfo = [];
		warning ('CDF_varInfo: Variable "%s" not found in "%s".', varName, fileInfo.Filename)
	end
end
