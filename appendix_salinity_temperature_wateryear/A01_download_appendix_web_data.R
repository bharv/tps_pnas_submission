# =============================================================================
# A01_download_appendix_web_data.R
# -----------------------------------------------------------------------------
# Pulls the web data used ONLY by the appendix:
#   * CDEC Delta Outflow Index, conductance and turbidity series   (Fig. S1, A08)
#   * CDEC + USGS NWIS conductance for the salinity contour maps    (Fig. S2, A09)
#   * TIGERweb hydrography (water mask for the contour maps)
# The functions live in 00_data/00_download_web_data.R; each skips files already
# on disk. Network time: ~10-20 min (USGS is requested month by month).
# =============================================================================

source("00_data/00_download_web_data.R")

download_tiger_hydro()
download_cdec_timeseries()       # -> data/web_downloads/cdec_timeseries/
download_cdec_salinity_field()   # -> data/web_downloads/salinity_field/cdec/
download_usgs_conductance()      # -> data/web_downloads/salinity_field/usgs_monthly/
