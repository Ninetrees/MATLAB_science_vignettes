% MEEdrift_writeCurrentPlotData

format long g

PlotLogFileID = fopen ([ SavePlotFilename, '.txt'], 'w');

fprintf (PlotLogFileID, '%s\n', mms_ql__EDI__BdvE__dataFile);
fprintf (PlotLogFileID, '%s\n', mms_ql__EDP__dataFile);

fprintf (PlotLogFileID, 'Beam record      : %d\n', iCenterBeam);
fprintf (PlotLogFileID, 'edi_gd_beam_t2k  : %d\n', centerBeam_t2k);
strCenterBeam = datestr (centerBeam_dn, 'yyyy-mm-dd HH:MM:ss.fff')
fprintf (PlotLogFileID, 'Beam time        : %s\n', strCenterBeam);

fprintf (PlotLogFileID, 'edp_E_interp     : %g %g %g\n', edp_E_interp);
fprintf (PlotLogFileID, 'B DMPA           : %g %g %g\n', centerBeamB);
fprintf (PlotLogFileID, 'B BPP            : %g %g %g\n', DMPA2BPP * centerBeamB);
fprintf (PlotLogFileID, '||B||            : %3.0f\n',    centerBeamB2n);

if S_star_dmpa (1) ~= -edi_d_dmpa_fillVal % S* is the negative of the drift step, so FILLVALs are negated, too
	fprintf (PlotLogFileID, 'S* BPP (converge): %g %g %g\n', S_star_bpp);
else
	fprintf (PlotLogFileID, 'S* BPP (converge): %g %g %g\n', S_star_dmpa);
end
fprintf (PlotLogFileID, 'Gyrofreq, period: %g Hz %g s\n', gyroFrequency, gyroPeriod);
fprintf (PlotLogFileID, 'DMPA2BPP\n');
fprintf (PlotLogFileID, '%16.13f %16.13f %16.13f\n', DMPA2BPP);

fclose (PlotLogFileID);
disp (['Beam analysis data saved as "', SavePlotFilename, '.txt', '".'])
