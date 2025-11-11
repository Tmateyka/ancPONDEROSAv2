import polars as pl
import pandas as pd
from pathlib import Path
from typing import Tuple, Union, List
import numpy as np

####### Load map file #######
class PlinkMap(pd.DataFrame):
    
    _metadata = ['map_len']
    
    def __init__(self, data: pd.DataFrame = None, *args, map_len: float = None, **kwargs):
        # Call parent DataFrame constructor with all args and kwargs
        if data is not None:
            # aDNA adaptation: normalize map on construction (sorted, unique bp, monotonic cm)
            data = self._sanitize_map_dataframe(data)
            super().__init__(data, *args, **kwargs)
        else:
            super().__init__(*args, **kwargs)
        
        # Only calculate genome length if not provided and we have the required columns
        if map_len is not None:
            self.map_len = map_len
        elif self._has_required_columns():
            self.map_len = self._calculate_genome_length()
        else:
            self.map_len = 0.0
    
    def _has_required_columns(self) -> bool:
        """Check if DataFrame has the required columns for genome length calculation."""
        required_cols = ['chromosome', 'cm']
        return all(col in self.columns for col in required_cols)

    @classmethod
    def from_file(cls, plink_file: Union[str, Path]):
        """Load PlinkMap from a plink file."""
        # aDNA adaptation: use wider integer types and tolerate extraneous whitespace/blank lines
        data = pd.read_csv(
            plink_file,
            sep=r'\s+',  # Any whitespace
            header=None,
            names=["chromosome", "rsid", "cm", "bp"],
            dtype={
                # int32 avoids potential overflow if chr codes are nontrivial in custom maps
                "chromosome": 'int32',
                "rsid": 'string',
                "cm": 'float64',
                "bp": 'int64'
            },
            comment='#',
            engine='python'
        )
        # aDNA adaptation: drop rows missing essential coords (rare but safer for sparse maps)
        data = data.dropna(subset=["chromosome", "cm", "bp"])
        return cls(data)

    def _calculate_genome_length(self) -> float:
        """Calculate total genome length from the map data."""
        if self.shape[0] == 0 or not self._has_required_columns():
            return 0.0
        
        tot = 0.0
        # Use regular pandas DataFrame operations to avoid recursion
        df = pd.DataFrame(self)  # Create a regular DataFrame copy
        # aDNA adaptation: ensure per-chrom sort before span calculation
        df = df.sort_values(["chromosome", "bp"], kind="mergesort")
        for _, chrom_df in df.groupby("chromosome", sort=False):
            if len(chrom_df) > 1:
                # ensure monotonic cm (already enforced in _sanitize_map_dataframe, but double-safe)
                cm = chrom_df["cm"].to_numpy()
                # span is last - first
                span = float(cm[-1] - cm[0])
                if span > 0:
                    tot += span
        return tot
    
    def get_col(self, col: str) -> np.ndarray:
        return self[col].to_numpy()

    @staticmethod
    def _sanitize_map_dataframe(df: pd.DataFrame) -> pd.DataFrame:
        """
        aDNA adaptation: make maps robust to sparse / noisy inputs.
        - sort by chromosome, bp
        - drop duplicate bp positions per chromosome (keep first)
        - enforce non-decreasing cm within each chromosome (cummax)
        """
        if df.empty:
            return df

        # Ensure expected columns exist; if not, return as-is
        required = {"chromosome", "bp", "cm"}
        if not required.issubset(df.columns):
            return df

        # Sort for stable processing
        df = df.sort_values(["chromosome", "bp"], kind="mergesort")

        # Drop duplicate bp per chromosome (keep first occurrence)
        df = df[~df.duplicated(subset=["chromosome", "bp"], keep="first")]

        # Enforce monotonic non-decreasing cm within each chromosome
        def _enforce_monotonic(group: pd.DataFrame) -> pd.DataFrame:
            cm = group["cm"].to_numpy()
            # Replace NaNs (if any) by forward fill, then cummax to prevent decreases
            # Forward fill
            if np.isnan(cm).any():
                # Simple forward fill for NumPy array
                mask = np.isnan(cm)
                # If first is NaN, set to 0 (conservative baseline)
                if mask[0]:
                    cm[0] = 0.0
                for i in range(1, cm.size):
                    if np.isnan(cm[i]):
                        cm[i] = cm[i-1]
            # Cummax to enforce non-decreasing
            np.maximum.accumulate(cm, out=cm)
            group = group.copy()
            group["cm"] = cm
            return group

        df = df.groupby("chromosome", group_keys=False, sort=False).apply(_enforce_monotonic)

        return df

    def filter_chrom(self, chrom: int) -> 'PlinkMap':
        # aDNA adaptation: return empty PlinkMap with correct columns if chrom missing
        filtered_df = self[self.chromosome == chrom]
        # Create new PlinkMap; map_len recalculated inside
        new_obj = PlinkMap(filtered_df)
        return new_obj    
    
    @property
    def _constructor(self):
        """Ensures that operations return PlinkMap objects instead of DataFrames."""
        def constructor(*args, **kwargs):
            # For internal pandas operations, pass map_len=0.0 to avoid recalculation
            # Extract map_len from kwargs if present, otherwise default to 0.0
            map_len = kwargs.pop('map_len', 0.0)
            return PlinkMap(*args, map_len=map_len, **kwargs)
        return constructor
    
    @property
    def _constructor_sliced(self):
        """Constructor for Series operations."""
        return pd.Series
    
    def __finalize__(self, other, method=None, **kwargs):
        """Propagate metadata during operations."""
        self = super().__finalize__(other, method, **kwargs)
        if hasattr(other, 'map_len'):
            self.map_len = other.map_len
        else:
            # For new objects created during operations, set map_len to 0
            # User can recalculate if needed
            self.map_len = 0.0
        return self
    
    @classmethod
    def from_dataframe(cls, df: pd.DataFrame, genome_length: float = None):
        """Create PlinkMap from existing DataFrame."""
        # aDNA adaptation: sanitize incoming frame to ensure stable interpolation
        df = cls._sanitize_map_dataframe(df)
        return cls(df, map_len=genome_length)
    
    def recalculate_map_length(self):
        """Explicitly recalculate the genome length."""
        self.map_len = self._calculate_genome_length()
        return self.map_len
    

class GeneticMap:

    def __init__(self, map_df: PlinkMap):

        self.map_df = map_df

        self.genome_len = map_df.map_len

    @classmethod
    def add_plink(cls, map_file: Union[str, Path]):

        map_df = PlinkMap.from_file(map_file)

        return cls(map_df)
    
    @classmethod
    def add_plink_list(cls, map_file_list: List[Union[str, Path]]):

        map_df_list = []
        for map_file in map_file_list:
            map_df = PlinkMap.from_file(map_file)
            map_df_list.append(map_df)

        map_df = PlinkMap.from_dataframe(pd.concat(map_df_list, ignore_index=True))

        return cls(map_df)

    @classmethod
    def add_hapmap(cls, hapmap_file: str):
        # Placeholder: not used in tests; left unchanged
        pass

    def _wget_hapmap(self, map_build: int) -> Tuple[pd.DataFrame, float]:
        # Placeholder: network fetch intentionally not implemented
        map_url = f"https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh{map_build}.map.zip"
        pass

    @classmethod
    def no_map(cls, map_build: int = 37):
        # Placeholder kept; no network I/O here
        # aDNA adaptation: For truly mapless runs, users should supply a plink map; we keep this stub to avoid surprises.
        raise NotImplementedError("no_map is not available offline. Provide a PLINK map file or per-chromosome map list.")

    def interp(self, arr: np.ndarray, chrom_arr: np.ndarray) -> np.ndarray:
        """
        Interpolate cM for given bp positions per chromosome.
        aDNA adaptation:
          - handle chromosomes with 0 or 1 map rows
          - tolerate unsorted/non-unique bp in map (sanitized earlier)
          - clamp outside range (np.interp behavior)
          - return NaN for chromosomes absent from map
        """
        arr = np.asarray(arr)
        chrom_arr = np.asarray(chrom_arr)

        interp_arr = np.full(arr.shape[0], np.nan, dtype=float)  # default NaN where impossible

        # Fast path: if map empty, return all NaN
        if self.map_df.shape[0] == 0:
            return interp_arr

        for chrom in np.unique(chrom_arr):
            chrom_idx = np.where(chrom_arr == chrom)[0]

            chrom_map = self.map_df.filter_chrom(chrom)

            # If no entries for this chromosome, leave NaNs
            if chrom_map.shape[0] == 0:
                continue

            cm_val = chrom_map.get_col("cm")
            bp_val = chrom_map.get_col("bp")

            # aDNA adaptation: ensure strictly increasing bp for interpolation
            # (sanitize again defensively)
            order = np.argsort(bp_val, kind="mergesort")
            bp_sorted = bp_val[order]
            cm_sorted = cm_val[order]

            # If only a single point exists, treat it as a flat map at that cM (clamp)
            if bp_sorted.size == 1:
                interp_arr[chrom_idx] = float(cm_sorted[0])
                continue

            # Drop any duplicated bp to satisfy np.interp precondition
            # (np.interp tolerates equal x? It requires ascending. We ensure unique.)
            uniq_mask = np.concatenate(([True], bp_sorted[1:] != bp_sorted[:-1]))
            bp_unique = bp_sorted[uniq_mask]
            cm_unique = cm_sorted[uniq_mask]

            if bp_unique.size == 1:
                # degenerate after de-duplication
                interp_arr[chrom_idx] = float(cm_unique[0])
                continue

            # Perform interpolation; np.interp clamps outside range to edge cm
            interp_arr[chrom_idx] = np.interp(arr[chrom_idx], bp_unique, cm_unique).astype(float)

        return interp_arr

    def get_genome_length(self):
        return self.genome_len
