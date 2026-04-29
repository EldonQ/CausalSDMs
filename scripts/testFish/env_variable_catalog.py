"""
env_variable_catalog.py
=======================
testFish | RiverATLAS v10

Scans ALL RiverATLAS v10 shapefiles to extract the complete list of
environmental variable columns, classifies them by ecological dimension,
and saves a reference catalog CSV for documentation and downstream use.

Output: output/RiverATLAS_v10_EnvVar_Catalog.csv
"""

from pathlib import Path
import calendar
import re

import geopandas as gpd
import pandas as pd

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent.resolve()
RIVER_DIR  = SCRIPT_DIR / "RiverATLAS_Data_v10_shp" / "RiverATLAS_v10_shp"
OUT_DIR    = SCRIPT_DIR / "output"
OUT_DIR.mkdir(exist_ok=True)

SHP_FILES = {
    "af":       RIVER_DIR / "RiverATLAS_v10_af.shp",
    "ar":       RIVER_DIR / "RiverATLAS_v10_ar.shp",
    "as":       RIVER_DIR / "RiverATLAS_v10_as.shp",
    "au":       RIVER_DIR / "RiverATLAS_v10_au.shp",
    "eu":       RIVER_DIR / "RiverATLAS_v10_eu.shp",
    "gr":       RIVER_DIR / "RiverATLAS_v10_gr.shp",
    "na":       RIVER_DIR / "RiverATLAS_v10_na.shp",
    "sa_north": RIVER_DIR / "RiverATLAS_v10_sa_north.shp",
    "sa_south": RIVER_DIR / "RiverATLAS_v10_sa_south.shp",
    "si":       RIVER_DIR / "RiverATLAS_v10_si.shp",
}

# Columns to exclude (not environmental variables)
META_EXCLUDE = {
    "geometry", "Hylak_id", "NextDownID", "MainRiver", "DistDown",
    "DistSea",  "CatchStorM", "VolR", "UpCumulS", "Shape_Leng", "Shape_Area",
    "geometry_type", "gid", "src_region",
}

# ── Category classification rules ───────────────────────────────────────────────
def classify(name: str) -> tuple[str, str]:
    """
    Classify a RiverATLAS column name into:
      (category, subcategory)

    Categories:
      Hydrology | Physiography | Climate | Landcover | Soils | Anthropogenic
    """
    n = name.lower()

    # Hydrology
    if any(k in n for k in ("dis_m3", "run_mm", "inu_pc", "lka_pc", "lkv_mc",
                              "rev_mc", "dor_pc", "ria_ha", "riv_tc", "gwt_cm")):
        return ("Hydrology", _hydrology_sub(n))

    # Physiography / Topography
    if any(k in n for k in ("ele_mt", "slp_dg", "sgr_dk")):
        return ("Physiography", "Elevation & Slope")

    # Climate
    if any(k in n for k in ("tmp_dc", "pre_mm", "pet_mm", "aet_mm",
                              "ari_ix", "cmi_ix", "snw_pc", "clz_cl", "cls_cl")):
        return ("Climate", _climate_sub(n))

    # Landcover
    if any(k in n for k in ("glc_cl", "glc_pc", "pnv_cl", "pnv_pc",
                              "wet_cl", "wet_pc", "for_pc", "crp_pc",
                              "pst_pc", "ire_pc", "gla_pc", "prm_pc",
                              "pac_pc", "tbi_cl", "tec_cl", "fmh_cl", "fec_cl")):
        return ("Landcover", _landcover_sub(n))

    # Soils & Geology
    if any(k in n for k in ("cly_pc", "slt_pc", "snd_pc", "soc_th",
                              "swc_pc", "lit_cl", "kar_pc", "ero_kh")):
        return ("Soils & Geology", _soils_sub(n))

    # Anthropogenic
    if any(k in n for k in ("pop_ct", "ppd_pk", "urb_pc", "nli_ix",
                              "rdd_mk", "hft_ix", "gad_id", "gdp_ud", "hdi_ix")):
        return ("Anthropogenic", _anthropogenic_sub(n))

    return ("Other", "Other")


def _hydrology_sub(n: str) -> str:
    if "dis_m3"  in n: return "Natural Discharge"
    if "run_mm"  in n: return "Land Surface Runoff"
    if "inu_pc"  in n: return "Inundation Extent"
    if "lka_pc"  in n: return "Lake Area (Limnicity)"
    if "lkv_mc"  in n: return "Lake Volume"
    if "rev_mc"  in n: return "Reservoir Volume"
    if "dor_pc"  in n: return "Degree of Regulation"
    if "ria_ha"  in n: return "River Area"
    if "riv_tc"  in n: return "River Volume"
    if "gwt_cm"  in n: return "Groundwater Table Depth"
    return "Hydrology"


def _climate_sub(n: str) -> str:
    if "tmp_dc"  in n: return "Air Temperature"
    if "pre_mm"  in n: return "Precipitation"
    if "pet_mm"  in n: return "Potential Evapotranspiration"
    if "aet_mm"  in n: return "Actual Evapotranspiration"
    if "ari_ix"  in n: return "Aridity Index"
    if "cmi_ix"  in n: return "Climate Moisture Index"
    if "snw_pc"  in n: return "Snow Cover"
    if "clz_cl"  in n: return "Climate Zone"
    if "cls_cl"  in n: return "Climate Stratum"
    return "Climate"


def _landcover_sub(n: str) -> str:
    if "glc_cl" in n or "glc_pc" in n: return "GLC2000 Land Cover"
    if "pnv_cl" in n or "pnv_pc" in n: return "Potential Natural Vegetation"
    if "wet_cl" in n or "wet_pc" in n: return "Wetlands"
    if "for_pc" in n: return "Forest Cover"
    if "crp_pc" in n: return "Cropland"
    if "pst_pc" in n: return "Pasture"
    if "ire_pc" in n: return "Irrigated Area"
    if "gla_pc" in n: return "Glacier"
    if "prm_pc" in n: return "Permafrost"
    if "pac_pc" in n: return "Protected Area"
    if "tbi_cl" in n: return "Terrestrial Biome"
    if "tec_cl" in n: return "Terrestrial Ecoregion"
    if "fmh_cl" in n: return "Freshwater Major Habitat"
    if "fec_cl" in n: return "Freshwater Ecoregion"
    return "Landcover"


def _soils_sub(n: str) -> str:
    if "cly_pc" in n: return "Clay Fraction"
    if "slt_pc" in n: return "Silt Fraction"
    if "snd_pc" in n: return "Sand Fraction"
    if "soc_th" in n: return "Soil Organic Carbon"
    if "swc_pc" in n: return "Soil Water Content"
    if "lit_cl" in n: return "Lithology"
    if "kar_pc" in n: return "Karst Area"
    if "ero_kh" in n: return "Soil Erosion"
    return "Soils & Geology"


def _anthropogenic_sub(n: str) -> str:
    if "pop_ct" in n: return "Population Count"
    if "ppd_pk" in n: return "Population Density"
    if "urb_pc" in n: return "Urban Extent"
    if "nli_ix" in n: return "Nighttime Lights"
    if "rdd_mk" in n: return "Road Density"
    if "hft_ix" in n: return "Human Footprint"
    if "gad_id" in n: return "GADM Admin Area"
    if "gdp_ud" in n: return "GDP"
    if "hdi_ix" in n: return "Human Development Index"
    return "Anthropogenic"


# ── Suffix explanations ────────────────────────────────────────────────────────
SUFFIX_EXPLAIN = {
    "_pyr":  "Annual average discharge at reach pour point",
    "_pmn":  "Annual minimum discharge at reach pour point",
    "_pmx":  "Annual maximum discharge at reach pour point",
    "_cyr":  "Annual average in reach catchment",
    "_cmn":  "Annual minimum in reach catchment",
    "_cmx":  "Annual maximum in reach catchment",
    "_cav":  "Average in reach catchment",
    "_cse":  "Spatial extent (%) in reach catchment",
    "_cmj":  "Spatial majority class in reach catchment",
    "_c01":  "January average in catchment",
    "_c02":  "February average in catchment",
    "_c03":  "March average in catchment",
    "_c04":  "April average in catchment",
    "_c05":  "May average in catchment",
    "_c06":  "June average in catchment",
    "_c07":  "July average in catchment",
    "_c08":  "August average in catchment",
    "_c09":  "September average in catchment",
    "_c10":  "October average in catchment",
    "_c11":  "November average in catchment",
    "_c12":  "December average in catchment",
    "_uav":  "Average in total upstream watershed",
    "_use":  "Spatial extent (%) in upstream watershed",
    "_usr":  "Sum in upstream watershed",
    "_ury":  "Annual average in upstream watershed",
    "_av":   "Average",
    "_su":   "Sum",
    "_mn":   "Minimum",
    "_mx":   "Maximum",
    "_yr":   "Annual",
    "_va":   "Value",
    "_lt":   "Long-term",
}


def explain_suffix(col: str) -> str:
    """Return a human-readable explanation of a column's suffix."""
    # Strip the prefix (first 6 chars, e.g. dis_m3)
    if len(col) <= 7:
        return ""
    suffix = col[7:]   # e.g. "_pyr", "_cyr", "_cmj"

    # Monthly
    m = re.match(r"^_c(\d{2})$", suffix)
    if m:
        month_num = int(m.group(1))
        month_name = calendar.month_abbr[month_num]
        return f"Monthly average in catchment ({month_name})"

    known = {
        "_pyr":  "Annual discharge at pour point",
        "_pmn":  "Annual min discharge at pour point",
        "_pmx":  "Annual max discharge at pour point",
        "_cyr":  "Annual average in catchment",
        "_cmn":  "Annual minimum in catchment",
        "_cmx":  "Annual maximum in catchment",
        "_cav":  "Average in catchment",
        "_cmj":  "Majority class in catchment",
        "_cse":  "Spatial extent in catchment (%)",
        "_uav":  "Average in upstream watershed",
        "_use":  "Spatial extent in upstream watershed (%)",
        "_usr":  "Sum in upstream watershed",
        "_ury":  "Annual average in upstream watershed",
        "_lt":   "Long-term",
        "_av":   "Average",
        "_su":   "Sum",
        "_mn":   "Minimum",
        "_mx":   "Maximum",
        "_yr":   "Annual",
        "_va":   "Value",
    }
    return known.get(suffix, f"Suffix: {suffix}")


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    print("=" * 65)
    print("  env_variable_catalog.py  |  RiverATLAS v10  |  Column Scanner")
    print("=" * 65)

    all_cols = set()
    first_shp = None

    print(f"\n[1] Scanning shapefiles ...")
    for region, shp_path in SHP_FILES.items():
        if not shp_path.exists():
            print(f"    [SKIP] {shp_path.name} — not found")
            continue
        print(f"    Scanning {shp_path.name} ...", end=" ", flush=True)
        try:
            gdf = gpd.read_file(str(shp_path), engine="pyogrio")
            print(f"{len(gdf):,} rows, {len(gdf.columns)} cols")
            for col in gdf.columns:
                if col not in META_EXCLUDE and col != "geometry":
                    all_cols.add(col)
            if first_shp is None:
                first_shp = gdf
        except Exception as e:
            print(f"ERROR: {e}")

    print(f"\n[2] Found {len(all_cols)} unique env/attribute columns")

    # ── Classify each column ─────────────────────────────────────────────────
    print("\n[3] Classifying columns by ecological dimension ...")

    records = []
    for col in sorted(all_cols):
        cat, subcat = classify(col)
        prefix = col.split("_")[0] if "_" in col else col[:6]
        suffix = explain_suffix(col)

        records.append({
            "column_name":     col,
            "prefix":         prefix,
            "category":       cat,
            "subcategory":    subcat,
            "suffix_explain": suffix,
            "dtypesample":    str(first_shp[col].dtype) if first_shp is not None else "N/A",
        })

    catalog_df = pd.DataFrame(records)

    # ── Summary ───────────────────────────────────────────────────────────────
    print("\n    Category summary:")
    for cat, grp in catalog_df.groupby("category"):
        print(f"      [{cat}] — {len(grp)} columns")

    # ── Save ─────────────────────────────────────────────────────────────────
    out_path = OUT_DIR / "RiverATLAS_v10_EnvVar_Catalog.csv"
    catalog_df.to_csv(out_path, index=False)
    print(f"\n[4] Saved: {out_path}")

    print("\n    Sample:")
    print(catalog_df[["column_name", "category", "subcategory"]].head(20).to_string(index=False))
    print("=" * 65)


if __name__ == "__main__":
    main()
