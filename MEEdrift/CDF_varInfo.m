function varInfo = CDF_varInfo (fileInfo, varName)

	% CDF_fileInfo_.Variables
	% Example (2015-07-24): 'Epoch' [1x2 double] [ 501] 'tt2000' 'T/'  'Full' 'None' [0] [-9223372036854775808]

	[nVariables, ~] = size (fileInfo.Variables);
	varFound = false;
	for iVariables = 1 : nVariables
		if strcmp (fileInfo.Variables {iVariables, 1}, varName)
			varInfo.name     = fileInfo.Variables {iVariables, 1};
			varInfo.dims     = fileInfo.Variables {iVariables, 2};
			varInfo.nRecs    = fileInfo.Variables {iVariables, 3};
			varInfo.type     = fileInfo.Variables {iVariables, 4};
			varInfo.dimsVary = fileInfo.Variables {iVariables, 5};
			varInfo.dummy1   = fileInfo.Variables {iVariables, 6};
			varInfo.dummy2   = fileInfo.Variables {iVariables, 7};
			varInfo.dummy3   = fileInfo.Variables {iVariables, 8};
			varInfo.padVal   = fileInfo.Variables {iVariables, 9};

			[nFILLVAL, ~] = size (fileInfo.VariableAttributes.FILLVAL);
			fillValFound = false;
			for iFILLVAL = 1 : nFILLVAL
				if strcmp (fileInfo.VariableAttributes.FILLVAL {iFILLVAL, 1}, varName)
					varInfo.fillVal  = fileInfo.VariableAttributes.FILLVAL {iFILLVAL, 2};
					fillValFound = true;
				end
			end
			if ~fillValFound
				varInfo.fillVal = [];
				warning ('CDF_varInfo: Variable "%s" FILLVAL not found in "%s".', varName, fileInfo.Filename)
			end

			[nSI_CONVERSION, ~] = size (fileInfo.VariableAttributes.SI_CONVERSION);
			SI_ConvsionFound = false;
			for iSI_CONVERSION = 1 : nSI_CONVERSION
				if strcmp (fileInfo.VariableAttributes.SI_CONVERSION {iSI_CONVERSION, 1}, varName)
					varInfo.SI_Convsion  = fileInfo.VariableAttributes.SI_CONVERSION {iSI_CONVERSION, 2};
					SI_ConvsionFound = true;
				end
			end
			if ~SI_ConvsionFound
				varInfo.SI_Convsion = [];
				warning ('CDF_varInfo: Variable "%s" SI_CONVERSION not found in "%s".', varName, fileInfo.Filename)
			end

			[nUNITS, ~] = size (fileInfo.VariableAttributes.UNITS);
			UnitsFound = false;
			for iUNITS = 1 : nUNITS
				if strcmp (fileInfo.VariableAttributes.UNITS {iUNITS, 1}, varName)
					varInfo.units  = fileInfo.VariableAttributes.UNITS {iUNITS, 2};
					UnitsFound = true;
				end
			end
			if ~UnitsFound
				varInfo.units = [];
				warning ('CDF_varInfo: Variable "%s" UNITS not found in "%s".', varName, fileInfo.Filename)
			end

			varFound = true;
		end
	end
	if ~varFound
		varInfo = [];
		warning ('CDF_varInfo: Variable "%s" not found in "%s".', varName, fileInfo.Filename)
	end
end
