% RecursiveMovingAverageFilter
function RMAF = RecursiveMovingAverageFilter (Data, FilterSize)
	% FilterSize should be odd > 1
	nData = uint32 (size (Data, 1)); % Expects a col of values, n x 1
	RMAF  = zeros (nData, 1);

	Accumulator = double (0.0);
	HalfWidth   = uint32 (FilterSize-1) / 2; % the point of this is to get (FilterSize - 1) / 2
  Accumulator = sum (Data (1: FilterSize)); % for FS = 3, this is 1:3

  RMAF (HalfWidth+1) = Accumulator; % for FS = 3, this is RMAF (2) = Acc

	i = uint32 (0);
	for i = HalfWidth+2: nData-HalfWidth  % for FS = 3, this is Acc = Acc - Data (3-1-1) + Data (4)
% 		if ~isnan (
% 		Accumulator = Accumulator - (Data (i-1-HalfWidth) * ~isnan (Data (i-1-HalfWidth))) + (Data (i+HalfWidth) * ~isnan (Data (i+HalfWidth)));
		Accumulator = Accumulator -  Data (i-1-HalfWidth)                                  +  Data (i+HalfWidth);
% 		if isnan (Accumulator)
% 			Accumulator = 0.0;
% 		end
		RMAF (i) = Accumulator;
	end
	RMAF = RMAF ./ FilterSize;

	RMAF (1: HalfWidth) = Data (1: HalfWidth);
	RMAF ((nData-HalfWidth+1): nData) = Data ((nData-HalfWidth+1): nData);
end
