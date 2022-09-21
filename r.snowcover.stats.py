#!/usr/bin/env python3

"""
MODULE:    r.snowcover.stats

AUTHOR(S): Stefan Blumentrath <stbl AT nve.no>

PURPOSE:   Compute maximum, mean and minimum snow cover from data provided by an ImageServer

COPYRIGHT: (C) 2022 NVE, Stefan Blumentrath, and by the GRASS Development Team

This program is free software under the GNU General Public License
(>=v2). Read the file COPYING that comes with GRASS for details.
"""

#%module
#% description: Compute maximum, mean and minimum snow cover from data provided by an ImageServer
#% keyword: raster
#% keyword: ImageServer
#% keyword: snow
#% keyword: Sentinel-3
#% keyword: DTM
#% keyword: zonal statistics
#%end

#%option
#% key: image_server
#% type: string
#% description: URL to ArcGIS ImageServer Typically starts with "http://"
#% answer: http://gis3.nve.no/image/rest/services/ImageService
#% guisection: Input
#% required: no
#%end

#%option
#% key: dtm_service
#% type: string
#% description: Name of the service of the ImageServer to query for DTM
#% answer: AuxDEM250
#% guisection: Input
#% multiple: no
#% required: no
#%end

#%option
#% key: snow_service
#% type: string
#% description: Name of the service of the ImageServer to query for snow cover data time series
#% answer: S3_SLSTR_fsc_sa
#% guisection: Input
#% multiple: no
#% required: no
#%end

#%option G_OPT_F_INPUT
#% key: aoi
#% type: string
#% description: Path to an OGR-readable file with Area of Interest
#% multiple: no
#% guisection: Input
#%end

#%option
#% key: date_start
#% type: string
#% description: First date (in ISO-format: YYYY-MM-DD) of the time series to compute snow cover for (defaults to yesterday)
#% required: no
#% guisection: Input
#%end

#%option
#% key: date_end
#% type: string
#% description: Last date (in ISO-format: YYYY-MM-DD) of the time series to compute snow cover for (defaults to today)
#% required: no
#% guisection: Input
#%end

#%option
#% key: snow_bands
#% type: integer
#% description: Altitude steps (bands) in meter to compute snow cover for
#% answer:250
#% guisection: Settings
#%end

#%option
#% key: min_snow_percent
#% type: integer
#% description: Minimum snow percent pre pixel to consider covered by snow
#% answer:20
#% guisection: Settings
#%end

# Currently not implemented
#%flag
#% key: d
#% description: Write temporary data to disk (default is to compute in memory) (comutation in memory is faster but may limit number of days that can be processed)
#% guisection: Settings
#%end

#%rules
#% requires: date_end, date_end
#%end

# snow_altitude_actinia aoi=C:\data\aoi.geojson
# Imports from builtin libs
# import atexit
from datetime import datetime, timedelta
import json
import logging
from math import ceil, floor
from pathlib import Path
from urllib import request, parse
import sys
import tempfile

# sys.path.insert(1, os.path.join(os.path.dirname(sys.path[0]), 'etc', 'r.in.wms'))

import grass.script as gscript

# Imports from non-standards Python libraries
# gdal må være tilgjengelig for Python interpreteren
# nedenfor defineres omgivelsesvariablene for CONDA

try:
    from osgeo import osr, gdal  # , ogr
except ImportError:
    gscript.error(
        _(
            "Cannot import gdal python library"
            "Please make sure that GDAL python bindings are installed (pip install GDAL)"
        )
    )
try:
    import numpy as np
except ImportError:
    gscript.error(
        _(
            "Cannot import numpy python library"
            "Please make sure that numpy installed (pip install numpy)"
        )
    )


def align_windows(window, ref):
    """Align two regions
    Python version of:
    https://github.com/OSGeo/grass/blob/main/lib/raster/align_window.c
    Modifies the input ``window`` to align to ``ref`` region. The
    resolutions in ``window`` are set to match those in ``ref``
    and the ``window`` edges (ymax, ymin, xmax, xmin) are modified
    to align with the grid of the ``ref`` region.
    The ``window`` may be enlarged if necessary to achieve the
    alignment. The ymax is rounded ymaxward, the ymin yminward,
    the xmax xmaxward and the xmin xminward. Lon-lon constraints are
    taken into consideration to make sure that the ymax doesn't go
    above 90 degrees (for lat/lon) or that the xmax does "wrap" past
    the xmin, etc.
    :param window: dict of window to align, with keys ymax, ymin, xmax,
                   xmin, pixelSizeY, pixelSizeX, is_latlong
    :type window: dict
    :param ref: dict of window to align to, with keys ymax, ymin, xmax,
                xmin, pixelSizeY, pixelSizeX, is_latlong
    :type ref: dict
    :return: a modified version of ``window`` that is aligend to ``ref``
    :rtype: dict
    """

    window["pixelSizeY"] = ref["pixelSizeY"]
    window["pixelSizeX"] = ref["pixelSizeX"]
    window["is_latlong"] = ref["is_latlong"]

    gscript.debug("before alignment:")
    gscript.debug(f'ymax: {window["ymax"]:.15g}')
    gscript.debug(f'ymin: {window["ymin"]:.15g}')
    gscript.debug(f'xmin: {window["xmin"]:.15g}')
    gscript.debug(f'xmax: {window["xmax"]:.15g}')

    window["ymax"] = (
        ref["ymax"]
        - floor((ref["ymax"] - window["ymax"]) / ref["pixelSizeY"]) * ref["pixelSizeY"]
    )
    window["ymin"] = (
        ref["ymin"]
        - ceil((ref["ymin"] - window["ymin"]) / ref["pixelSizeY"]) * ref["pixelSizeY"]
    )
    # Rast_easting_to_col() wraps easting:
    # xmax can become < xmin, or both xmin and xmax are shifted */
    window["xmin"] = (
        ref["xmin"]
        + floor((window["xmin"] - ref["xmin"]) / ref["pixelSizeX"]) * ref["pixelSizeX"]
    )
    window["xmax"] = (
        ref["xmax"]
        + ceil((window["xmax"] - ref["xmax"]) / ref["pixelSizeX"]) * ref["pixelSizeX"]
    )

    if window["is_latlong"]:
        while window["ymax"] > 90.0 + window["pixelSizeY"] / 2.0:
            window["ymax"] -= window["pixelSizeY"]
        while window["ymin"] < -90.0 - window["pixelSizeY"] / 2.0:
            window["ymin"] += window["pixelSizeY"]

    window["width"] = int((window["xmax"] - window["xmin"]) / window["pixelSizeX"])
    window["height"] = int((window["ymax"] - window["ymin"]) / window["pixelSizeY"])
    window[
        "bbox"
    ] = f"{window['xmin']}, {window['ymin']}, {window['xmax']}, {window['ymax']}"

    gscript.debug("after alignment:")
    gscript.debug(f'ymax: {window["ymax"]:.15g}')
    gscript.debug(f'ymin: {window["ymin"]:.15g}')
    gscript.debug(f'xmin: {window["xmin"]:.15g}')
    gscript.debug(f'xmax: {window["xmax"]:.15g}')

    return window


def create_gdal_raster(
    reference_grid, srs, data_type, ds_name="dtm_raster", driver_name="MEM"
):
    """Create an empty GDAL raster dataset to work with"""
    driver = gdal.GetDriverByName(driver_name)
    target_ds = driver.Create(
        ds_name,
        int(
            (reference_grid["xmax"] - reference_grid["xmin"])
            / reference_grid["pixelSizeX"]
        ),
        int(
            (reference_grid["ymax"] - reference_grid["ymin"])
            / reference_grid["pixelSizeY"]
        ),
        1,
        data_type,
    )
    target_ds.GetRasterBand(1).Fill(0)
    target_ds.SetGeoTransform(
        (
            (reference_grid["xmin"]),
            reference_grid["pixelSizeX"],
            0,
            reference_grid["ymax"],
            0,
            -reference_grid["pixelSizeY"],
        )
    )
    target_ds.SetProjection(srs.ExportToWkt())
    return target_ds


def aggregate_raster(
    in_ds,
    reference_grid,
    spatial_ref,
    data_type=gdal.GDT_Float32,
    agg_alg=gdal.GRA_Average,
):
    """Aggregate Raster to coarser resolution"""
    target_ds = create_gdal_raster(reference_grid, spatial_ref, data_type)
    gdal.Warp(
        target_ds,
        in_ds,
        # format=gdal_format,
        # width=ref_grid['width'],
        # height=ref_grid['height'],
        resampleAlg=agg_alg,
        # dstSRS=spatial_ref,
        multithread=True,
    )

    out_array = np.array(target_ds.ReadAsArray())
    target_ds = None
    return out_array


def reproject_geojson(geo_json, spatial_ref):
    """Reproject GeoJSON to CRS of Image server"""
    # Set GeoJSON spatial reference
    s_srs = osr.SpatialReference()
    s_srs.ImportFromEPSG(4326)

    # Create temporary (shape) file
    tmp_shp = tempfile.NamedTemporaryFile().name + ".shp"

    # Reproject AOI layer
    gdal.VectorTranslate(
        tmp_shp, geo_json, dstSRS=spatial_ref, srcSRS=s_srs, reproject=True
    )

    return gdal.OpenEx(tmp_shp)


def rasterize_aoi(reference_grid, aoi, spatial_ref):
    """Rasterize AOI on reference grid"""

    # Create target raster map from reference grid
    target_ds = create_gdal_raster(
        reference_grid,
        spatial_ref,
        gdal.GDT_Byte,
        ds_name="mask_raster",
        driver_name="MEM",
    )

    # Rasterize AOI on target raster map
    gdal.Rasterize(
        target_ds, aoi, allTouched=True, layers=[aoi.GetLayerByIndex(0).GetName()]
    )
    return np.array(target_ds.ReadAsArray(), dtype=np.bool_)


def query_image_server(
    image_service="http://gis3.nve.no/image/rest/services/ImageService/S3_SLSTR_fsc_sa/ImageServer/",
    query_params=None,
    print_metadata=None,
    aoi=None,
):
    """Query an ImageServer instance for data in space and time"""
    # ToDo: Probably more elegantly implemented as a class
    # (checking server, listing services, ...)
    # Get Service metadata
    url = image_service + "?f=pjson"
    with request.urlopen(url) as response:
        if not response.getcode() != "200":
            gscript.error(response.read())
            sys.exit("Kan ikke lese fra ImageServeren.")
        response_text = response.read()
    service_description = json.loads(response_text.decode("utf-8"))
    if print_metadata and print_metadata == "service_description":
        return service_description

    # Get IDs of image(s) for selected day(s)
    params = {
        "returnGeometry": "false",
        "returnIdsOnly": "true",
        "f": "json",
    }
    if query_params:
        params.update(query_params)

    query_string = parse.urlencode(params)

    url = image_service + "query?" + query_string

    with request.urlopen(url) as response:
        response_text = response.read()
    j = json.loads(response_text.decode("utf-8"))
    if "objectIds" not in j:
        logging.error("Feil med spørring til ImageServeren.")
    images = j["objectIds"]
    if not images:
        gscript.warning("Ingen bilder på ImageServeren for valgt tidsrom.")

    # Get url of image(s) for selected day
    if not aoi:
        aoi = service_description["fullExtent"]
        aoi["pixelSizeX"] = service_description["pixelSizeX"]
        aoi["pixelSizeY"] = service_description["pixelSizeY"]

    width, height = (
        int((aoi["xmax"] - aoi["xmin"]) / aoi["pixelSizeX"]),
        int((aoi["ymax"] - aoi["ymin"]) / aoi["pixelSizeY"]),
    )
    bbox = ",".join(map(str, [aoi["xmin"], aoi["ymin"], aoi["xmax"], aoi["ymax"]]))
    image_params = {
        "bbox": bbox,
        "format": "tiff",
        "pixelType": service_description["pixelType"],
        "noData": "",
        "size": f"{width},{height}",
        "interpolation": "RSP_NearestNeighbor",
        "f": "json",
    }
    query_string = parse.urlencode(image_params)

    image_links = {}
    for i in images:
        url = image_service + str(i) + "/image?" + query_string
        with request.urlopen(url) as response:
            response_text = response.read()
        j = json.loads(response_text.decode("utf-8"))
        image_links[i] = j["href"]

    return image_links


def get_snow_statistics_for_height_interval(
    snow_min_elevation, snow_max_elevation, np_dtm, np_snow_ma, class_nr
):
    """Compute snow statistics for altitude range"""
    gscript.verbose("Beregner verdier for klasse {}".format(class_nr))
    gscript.verbose("Klassens minimumshøyde er {}m".format(snow_min_elevation))
    gscript.verbose("Klassens maximumshøyde er {}m".format(snow_max_elevation))

    np_dtm_class_ma = np.ma.masked_where(
        (np_dtm < snow_min_elevation) | (np_dtm > snow_max_elevation), np_dtm, copy=True
    )

    mask = np.ma.getmask(np_dtm_class_ma)
    np_snow_ma_ma = np.ma.MaskedArray(np_snow_ma, mask)

    snow_class_mean_percentage = round(np_snow_ma_ma.mean(), 2) or -9999
    snow_class_min_percentage = round(np_snow_ma_ma.min(), 2) or -9999
    snow_class_max_percentage = round(np_snow_ma_ma.max(), 2) or -9999

    gscript.verbose(
        "Klassens gjennomsnittlig snøprosent er {}%".format(snow_class_mean_percentage)
    )
    gscript.verbose(
        "Klassens minimum snøprosent er {}%".format(snow_class_min_percentage)
    )
    gscript.verbose(
        "Klassens maximum snøprosent er {}%".format(snow_class_max_percentage)
    )

    return (
        snow_class_mean_percentage,
        snow_class_min_percentage,
        snow_class_max_percentage,
    )


DATA_TYPES = {
    "U1": gdal.GDT_Byte,
    "U2": gdal.GDT_Byte,
    "U4": gdal.GDT_Byte,
    "U8": gdal.GDT_Byte,
    "U16": gdal.GDT_UInt16,
    "U32": gdal.GDT_UInt32,
    "S8": gdal.GDT_Byte,
    "S16": gdal.GDT_Int16,
    "S32": gdal.GDT_Int32,
    "F32": gdal.GDT_Float32,
    "F64": gdal.GDT_Float64,
    "C64": gdal.GDT_CFloat64,
    "C128": gdal.GDT_CFloat64,
    # gdal.GDT_CFloat32
    # gdal.GDT_CFloat64
    # gdal.GDT_CInt16
    # gdal.GDT_CInt32
    "UNKNOWN": gdal.GDT_Unknown,
}


def main():
    """Do the main work"""
    # User defined variables
    aoi = options["aoi"]
    server = options["image_server"]
    dtm_service = options["dtm_service"]
    snow_service = options["snow_service"]
    date_start = options["date_start"]
    date_end = options["date_start"]

    if not date_start:
        date_start = (datetime.now() - timedelta(days=1)).strftime("%Y-%m-%d")

    min_snow_percent = int(options["min_snow_percent"])
    # Probably more userfriendly for generate relatable height intervals dynamically and
    # let user define height intervals here
    snow_bands = int(options["snow_bands"])

    # Currently not used
    # not_in_memory = flags["d"]

    if any([protocol in aoi for protocol in ["http", "ftp"]]):
        aoi = f"/vsicurl/{aoi}"
    elif not Path(aoi).exists():
        gscript.fatal(_("Finner ikke {}. Sjekk filnavn og sti.".format(aoi)))

    try:
        dt_date_start = datetime.fromisoformat(date_start)
    except ValueError:
        gscript.fatal(
            _("Parameter 'start dato' er ikke et ISO-formatert dato (YYYY-MM-DD).")
        )

    if not date_end:
        dt_date_end = datetime.today()
        date_end = dt_date_end.strftime("%Y-%m-%d")
    else:
        try:
            dt_date_end = datetime.fromisoformat(date_end)
        except ValueError:
            gscript.fatal(
                _("Parameter 'end dato' er ikke et ISO-formatert dato (YYYY-MM-DD).")
            )
    if dt_date_start > dt_date_end:
        gscript.fatal(_("Valgt start dato ligger etter valgt slutt dato."))

    # Compile URL and query paramters for snow image server
    query = {
        "where": f"""opptakstidspunkt >= date '{date_start} 00:00:00' AND opptakstidspunkt <= date '{date_end} 23:59:59'"""
    }
    gscript.debug("SQL er: {query['where']}")

    image_server = "{server}/{service}/ImageServer/"

    snow_image_service = image_server.format(server=server, service=snow_service)
    snow_service_description = query_image_server(
        image_service=snow_image_service, print_metadata="service_description"
    )

    # Get spatial reference of ImageServer
    spatial_ref = osr.SpatialReference()
    spatial_ref.ImportFromEPSG(
        snow_service_description["fullExtent"]["spatialReference"]["wkid"]
    )

    # Get Area of interest
    aoi_reproj = reproject_geojson(aoi, spatial_ref)
    aoi_layer = aoi_reproj.GetLayerByIndex(0)
    aoi_dict = {}
    (
        aoi_dict["xmin"],
        aoi_dict["xmax"],
        aoi_dict["ymin"],
        aoi_dict["ymax"],
        aoi_dict["is_latlong"],
    ) = (*aoi_layer.GetExtent(), bool(spatial_ref.IsGeographic()))

    # Create a reference grid to operate on (extent, resolution, rows, cols)
    ref_grid = snow_service_description["fullExtent"]
    ref_grid["is_latlong"] = bool(spatial_ref.IsGeographic())
    ref_grid["pixelSizeX"], ref_grid["pixelSizeY"] = (
        snow_service_description["pixelSizeX"],
        snow_service_description["pixelSizeY"],
    )

    # Align Area of Interest to reference grid
    raster_aoi = align_windows(aoi_dict, ref_grid)

    # Rasterize Area of interest on reference grid to be sed as mask
    raster_mask = rasterize_aoi(raster_aoi, aoi_reproj, spatial_ref)

    # Get a list of images for the requested time period from the server
    images = query_image_server(
        image_service=snow_image_service, query_params=query, aoi=raster_aoi
    )

    # Create an empty numpy array
    np_snow = np.zeros(raster_mask.shape)
    np_snow_valid_count = np.zeros(raster_mask.shape, dtype=np.byte)

    # Read images to array and compute average value over available images
    for i in images.values():
        img_ds = gdal.Open(f"/vsicurl/{i}")
        img_ds = np.array(img_ds.ReadAsArray())
        np_snow += np.where((img_ds > 100) & (img_ds <= 200), ds - 100, 0)
        np_snow_valid_count += np.where((img_ds > 100) & (img_ds <= 200), 1, 0).astype(np.bool_)
        img_ds = None
    np_snow = np.divide(np_snow, np_snow_valid_count, where=np_snow_valid_count != 0)

    try:
        # Ekstra test siden extract by mask ikke feiler i Desktop
        # bare pixler fra og med brukerdefinert grenseverdi skal være med
        np_snow_ma = np.ma.masked_where(np_snow < min_snow_percent, np_snow, copy=True)
        if np_snow_ma.count() == 0:
            gscript.fatal(
                _(
                    "Fant ingen gyldige piksler. Dette er enten fordi det ikke finnes raster for valgt dato eller at angitt"
                    "polygon ikke overlapper med raster."
                )
            )

        px_area = raster_aoi["pixelSizeX"] * raster_aoi["pixelSizeY"]

        res_dict = {
            "snow_area": float(np_snow_ma.count() * px_area / 1000000),
            "snow_mean_percentage": float(round(np_snow_ma.mean(), 2)),
            "snow_min_percentage": float(round(np_snow_ma.min(), 2)),
            "snow_max_percentage": float(round(np_snow_ma.max(), 2)),
        }

        gscript.verbose(
            "Gjennomsnittlig snøprosent er {}%".format(res_dict["snow_mean_percentage"])
        )
        gscript.verbose(
            "Minimum snøprosent er {}%".format(res_dict["snow_min_percentage"])
        )
        gscript.verbose(
            "Maximum snøprosent er {}%".format(res_dict["snow_max_percentage"])
        )

        # np_snow_ma_non_zero = np.ma.masked_where(np_snow > 254, np_snow, copy=True)
        res_dict["area_complete"] = float(np_snow_ma.count() * px_area / 1000000)

        gscript.verbose(
            "Areal for piklser med minst {min_snow}% snø er {areal}km2".format(
                min_snow=min_snow_percent, areal=res_dict["snow_area"]
            )
        )
        gscript.verbose("Analyseareal er {}km2".format(res_dict["area_complete"]))

        res_dict["area_with_snow"] = float(
            round(res_dict["snow_area"] / res_dict["area_complete"] * 100, 2)
        )

        # Find elevation for pixels with snow
        res_dict["snow_mean_elevation"] = -9999
        res_dict["snow_min_elevation"] = -9999
        res_dict["snow_max_elevation"] = -9999

        # Compute on altitude bands
        try:
            if (res_dict["snow_min_percentage"]) > 0:

                dtm_image_service = image_server.format(
                    server=server, service=dtm_service
                )
                dtm_service_description = query_image_server(
                    image_service=dtm_image_service,
                    print_metadata="service_description",
                )
                dtm_ref_grid = dtm_service_description["fullExtent"]
                dtm_ref_grid["is_latlong"] = bool(spatial_ref.IsGeographic())
                dtm_ref_grid["pixelSizeX"], dtm_ref_grid["pixelSizeY"] = (
                    dtm_service_description["pixelSizeX"],
                    dtm_service_description["pixelSizeY"],
                )
                dtm_aoi = align_windows(aoi_dict.copy(), dtm_ref_grid)

                dtm = query_image_server(image_service=dtm_image_service, aoi=dtm_aoi)
                np_dtm = aggregate_raster(
                    gdal.Open(f"/vsicurl/{list(dtm.values())[0]}"),
                    raster_aoi,
                    spatial_ref,
                    agg_alg=gdal.GRA_Average,
                )

                # Create snow mask
                mask = np.ma.getmask(np_snow_ma)
                # Mask DTM
                np_dtm_ma = np.ma.MaskedArray(np_dtm, mask)

                res_dict["snow_mean_elevation"] = float(round(np_dtm_ma.mean(), 2))
                res_dict["snow_min_elevation"] = float(round(np_dtm_ma.min(), 2))
                res_dict["snow_max_elevation"] = float(round(np_dtm_ma.max(), 2))

                gscript.verbose(
                    "Gjennomsnittlig terrenghøyde med snø er {} meter".format(
                        res_dict["snow_mean_elevation"]
                    )
                )
                gscript.verbose(
                    "Laveste terrenghøyde med snø er {} meter".format(
                        res_dict["snow_min_elevation"]
                    )
                )
                gscript.verbose(
                    "Høyeste terrenghøyde med snø er {} meter".format(
                        res_dict["snow_max_elevation"]
                    )
                )

                snow_classes = [
                    snow_class
                    for snow_class in np.unique(
                        np.floor(np_dtm_ma / snow_bands).astype(np.int8)
                    )
                    if snow_class
                ]

                snow_diff_elevation = (
                    res_dict["snow_max_elevation"] - res_dict["snow_min_elevation"]
                ) / len(snow_classes)
                gscript.verbose(
                    "snow_diff_elevation er {} meter".format(snow_diff_elevation)
                )

                for i in snow_classes:
                    snow_class_min = i * snow_bands
                    snow_class_max = snow_class_min + snow_bands
                    (
                        snow_class_mean_percentage,
                        snow_class_min_percentage,
                        snow_class_max_percentage,
                    ) = get_snow_statistics_for_height_interval(
                        snow_class_min, snow_class_max, np_dtm, np_snow_ma, i
                    )

                    res_dict[
                        f"snow_{snow_class_min}m_to_{snow_class_max}m_mean_percentage"
                    ] = float(snow_class_mean_percentage)

                    res_dict[
                        f"snow_{snow_class_min}m_to_{snow_class_max}m_min_percentage"
                    ] = float(snow_class_min_percentage)

                    res_dict[
                        f"snow_{snow_class_min}m_to_{snow_class_max}m_max_percentage"
                    ] = float(snow_class_max_percentage)
            else:
                gscript.warning(
                    "Minimum snøprosent er {}. Beregner derfor ikke høyde".format(
                        res_dict["snow_min_percentage"]
                    )
                )
        except RuntimeError as err:
            gscript.fatal("Feilet ved DTM-analyse: {}".format(err))

        print(json.dumps(res_dict))
        # return json.dumps(res_dict)

    except RuntimeError as err:
        gscript.fatal("Feilet: {}".format(err))


if __name__ == "__main__":
    options, flags = gscript.parser()
    # atexit.register(cleanup)
    sys.exit(main())
