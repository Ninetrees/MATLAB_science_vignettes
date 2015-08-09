% MEEdrift_plot_EDP_zoomed_region
% disp 'MEEdrift_plot_EDP_zoomed_region' % V&V

if ~exist ('hEDP_subplot_zoomAxes')
	hEDP_subplot_zoomAxes = subplot (1, 8, [6,8]);
	EDP_plotZoomPos = get (hEDP_subplot_zoomAxes, 'position'); % [left,bottom,width,height]
	set (hEDP_subplot_zoomAxes, 'position', ...
		[EDP_plotZoomPos(1)+0.03 EDP_plotZoomPos(2) EDP_plotZoomPos(3) EDP_plotZoomPos(4)]);
else
	axes (hEDP_subplot_zoomAxes);
	delete (hEDP_plotZoom);
end

zoomTimeSeries {1} = edp_dn (edp_zoomStart: edp_zoomEnd);
zoomDataSeries {1} = edp_E3D_dsl (:, edp_zoomStart: edp_zoomEnd);
if FilterData
	zoomDataSeriesf {1} = edp_E3D_dslf (:, edp_zoomStart: edp_zoomEnd);
else
	zoomDataSeriesf {1} = edp_E3D_dsl (:, edp_zoomStart: edp_zoomEnd);
end

hEDP_plotZoom = plot (hEDP_subplot_zoomAxes, ...
	zoomTimeSeries {1}, zoomDataSeries {1}, ...
	zoomTimeSeries {1}, zoomDataSeriesf {1}, ...
	'LineStyle', '-', 'Marker', '.', 'MarkerSize', 2);

zoom_dateStart = timeSeries {1} (edp_zoomStart);
zoom_dateEnd   = timeSeries {1} (edp_zoomEnd);
xlim ([zoom_dateStart zoom_dateEnd]);
xTickData = linspace (zoom_dateStart, zoom_dateEnd, 5);
set (gca, 'XTick', xTickData);
datetick ('x', 'MM:SS','keepticks', 'keeplimits');
ZoomDateTicks ('on');
