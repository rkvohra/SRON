# SRON
Scripts made during internship a SRON which analyse correlation between agricutlural burning and Carbon Monoxide (CO) levels in Northern India.
"create..." scripts import data (which is not available due to confidentiality) from different sources:
 - Fires data from VIIRS, GFED, GFAS, SAGE
 - WRF model simulation data of CO in the region
And process it and export it as JSON.
"timeseries.py" further processes data based on location of interest and plots time series with some statistical measures over desired time period.
"barplot.py" does the same but for bar plots of cumulative burning.
