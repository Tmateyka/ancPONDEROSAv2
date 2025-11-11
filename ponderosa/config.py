"""
Configuration management for PONDEROSA genetic relationship inference.

This module provides dataclass-based configuration management for file inputs,
algorithm parameters, and output settings. Supports loading from YAML files
and command-line arguments.

Classes:
    FilesConfig: Configuration for input files
    AlgorithmConfig: Configuration for algorithm parameters
    OutputConfig: Configuration for output settings
    PonderosaConfig: Main configuration container
    
Example:
    >>> config = PonderosaConfig.from_yaml("ponderosa.yaml")
    >>> config.validate()
    >>> analyzer = PonderosaAnalyzer(config)
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Union, Dict, Any
import yaml
import logging

__all__ = ['FilesConfig', 'AlgorithmConfig', 'OutputConfig', 'PonderosaConfig']

logger = logging.getLogger(__name__)


@dataclass
class FilesConfig:
    """File input configuration for PONDEROSA."""
    
    # Required files
    fam: Path
    ibd_caller: str
    
    # Optional files
    ages: Optional[Path] = None
    populations: Optional[Path] = None
    training: Optional[Path] = None
    rel_tree: Optional[Path] = None
    priors: Optional[Path] = None

    mapf: Optional[Path] = None
    map_files: Optional[List[Path]] = None
    _map_file_list: Optional[List[Path]] = field(default=None, init=False)

    ibd: Optional[Path] = None
    ibd_files: Optional[List[Path]] = None
    _ibd_file_list: Optional[List[Path]] = field(default=None, init=False)

    # aDNA: optional inputs for low-coverage/inbred ancient samples
    kin_file: Optional[Path] = None             # e.g., lcMLkin/KIN pairwise output
    read_file: Optional[Path] = None            # READ pairwise output
    haproh_file: Optional[Path] = None          # hapROH per-individual summary
    haproh_files: Optional[List[Path]] = None   # allow multiple hapROH files

    def _validate_file_lists(self, single_file: Path, file_list: List[Path], file_type: str) -> List[Path]:

        single_file_arg = lambda x: {"map": "mapf", "ibd": "ibd"}[x]
        multiple_file_arg = lambda x: f"{x}_files"

        if single_file and file_list:
            raise ValueError(f"Cannot specify both '{single_file_arg(file_type)}' and '{multiple_file_arg(file_type)}'. Choose one.")

        if single_file:
            if not single_file.exists():
                raise FileNotFoundError(f"{single_file_arg(file_type)} file not found: {single_file}")
            
            if "chr1" in single_file.name: # Populate with all other chromosomes
                out_list = [single_file.with_name(single_file.name.replace("chr1", f"chr{i}")) for i in range(1, 23)]
            else:
                out_list = [single_file]

        elif file_list:
            if not file_list:
                raise ValueError(f"{multiple_file_arg(file_type)} list cannot be empty")
            out_list = []
            for i, file_n in enumerate(file_list):
                if not file_n.exists():
                    raise FileNotFoundError(f"{multiple_file_arg(file_type)} file not found: {file_n} (index {i})")
                out_list.append(file_n)

        else:
            # aDNA: in KIN/READ-only mode there may be no IBD or MAP files; caller controls when we use this helper.
            # leaving behavior unchanged for the legacy segments pipeline (we do NOT raise here to avoid changing flow)
            raise ValueError(f"{file_type} file(s) have not been provided. Specify either '{single_file_arg(file_type)}' or '{multiple_file_arg(file_type)}'.")  # FIX: ensure this actually raises

        return out_list
    
    def validate(self) -> None:
        """Validate that required files exist."""
        # Check required files
        required_files = [self.fam]
        for file_path in required_files:
            if not file_path.exists():
                raise FileNotFoundError(f"Required file not found: {file_path}")
        
        # Determine input mode(s)
        has_ibd = bool(self.ibd or self.ibd_files)
        has_kin_read = bool(self.kin_file or self.read_file)  # aDNA

        # aDNA: validate KIN/READ/hapROH paths if provided (non-fatal warnings if missing)
        for optional_path, label in [
            (self.kin_file, "kin_file"),
            (self.read_file, "read_file"),
            (self.haproh_file, "haproh_file"),
        ]:
            if optional_path and not Path(optional_path).exists():
                logger.warning(f"Optional aDNA file not found: {label} -> {optional_path}")

        if self.haproh_files:
            for i, hp in enumerate(self.haproh_files):
                if not Path(hp).exists():
                    logger.warning(f"Optional aDNA file not found: haproh_files[{i}] -> {hp}")
        
        # Legacy segments-based pipeline (unchanged): require IBD + MAP (MAP strongly recommended for hap-ibd)
        if has_ibd:
            self._ibd_file_list = self._validate_file_lists(self.ibd, self.ibd_files, "ibd")
            # Only require maps when a map-driven caller is used; we keep parity with your tests which supply mapf
            try:
                self._map_file_list = self._validate_file_lists(self.mapf, self.map_files, "map")
            except ValueError as e:
                # If the caller is hap-ibd, maps are expected; otherwise allow missing maps.
                caller = (self.ibd_caller or "").lower()
                if caller in {"hap-ibd", "hapibd"}:
                    raise
                else:
                    # no map needed for callers that emit cM lengths directly
                    self._map_file_list = []
        else:
            # aDNA: allow KIN/READ-only mode (no IBD segments), typically used for ultra-low coverage/pseudo-haploid data
            if has_kin_read:
                self._ibd_file_list = []
                self._map_file_list = []
                logger.info("aDNA mode: proceeding without IBD segment files (KIN/READ provided).")
            else:
                # No IBD and no aDNA surrogates => nothing to analyze
                raise ValueError(
                    "No IBD files provided and no aDNA surrogates (kin/read). "
                    "Provide 'ibd'/'ibd_files' OR 'kin_file'/'read_file'."
                )

        # Check optional files if provided
        optional_files = [self.ages, self.mapf, self.populations, self.training, self.priors, self.rel_tree]
        for file_path in optional_files:
            if file_path and not Path(file_path).exists():
                logger.warning(f"Optional file not found: {file_path}")

    @property
    def map_file_list(self) -> List[Path]:
        """Return unified list of map files. Must call validate() first."""
        if self._map_file_list is None:
            raise RuntimeError("Must call validate() before accessing map_file_list")
        return self._map_file_list
    
    @property
    def ibd_file_list(self) -> List[Path]:
        """Return unified list of map files. Must call validate() first."""
        if self._ibd_file_list is None:
            raise RuntimeError("Must call validate() before accessing map_file_list")
        return self._ibd_file_list

@dataclass 
class AlgorithmConfig:
    """Algorithm parameter configuration for PONDEROSA."""
    
    # IBD processing parameters
    min_segment_length: float = 3.0
    min_total_ibd: float = 50.0
    max_gap: float = 1.0
    
    # Relationship inference parameters
    use_phase_correction: bool = True
    mean_phase_error_distance: float = 50.0
    
    # Classification thresholds
    degree_threshold: float = 0.8
    haplotype_threshold: float = 0.2
    
    # Population parameters
    population: str = "pop1"
    genome_length: float = 3545.0  # Default human genome length in cM
    
    # Performance parameters
    validate_pair_order: bool = True
    parallel_processing: bool = True

    # aDNA: knobs for degraded/low-coverage data
    pseudo_haploid: bool = False                 # treat data as pseudo-haploid if caller prepared it that way
    min_overlap_snps: int = 0                    # minimum shared SNPs for pair to be considered (for KIN/READ)
    trim_damage_bases: Optional[int] = None      # soft-trim deaminated ends if upstream didn’t
    min_pairs_for_normalization: int = 0         # READ normalization cohort size (0 = no extra requirement)
    allow_no_segments: bool = True               # permit KIN/READ-only workflows

@dataclass
class OutputConfig:
    """Output configuration for PONDEROSA."""
    
    # Output settings
    output: str = "ponderosa_results"
    min_probability: float = 0.5
    
    # Output formats
    write_readable: bool = True
    write_pickle: bool = True
    write_detailed: bool = False
    write_training: bool = False
    
    # Visualization options
    create_plots: bool = False
    plot_format: str = "png"
    plot_dpi: int = 300
    
    # Logging
    log_level: str = "INFO"
    verbose: bool = False

@dataclass
class PonderosaConfig:
    """Main PONDEROSA configuration container."""
    
    files: FilesConfig
    algorithm: AlgorithmConfig = field(default_factory=AlgorithmConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    
    @classmethod
    def from_yaml(cls, yaml_path: Union[str, Path]) -> "PonderosaConfig":
        """Load configuration from YAML file."""
        yaml_path = Path(yaml_path)
        
        if not yaml_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {yaml_path}")
        
        with open(yaml_path) as f:
            config_dict = yaml.safe_load(f)
        
        return cls.from_dict(config_dict)
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> "PonderosaConfig":
        """Create configuration from dictionary (CLI args or YAML)."""
        
        files_dict = config_dict.get("files", {})
        
        # Single path fields (convert to Path objects)
        single_path_fields = {
            'ibd', 'fam', 'ages', 'mapf', 'populations', 'training', 'rel_tree',
            # aDNA:
            'kin_file', 'read_file', 'haproh_file'
        }
        for key, value in list(files_dict.items()):
            if value is not None and key in single_path_fields:
                files_dict[key] = Path(value)
        
        # Handle list fields (convert each string to Path)
        if 'map_files' in files_dict and files_dict['map_files'] is not None:
            files_dict['map_files'] = [Path(f) for f in files_dict['map_files']]
        if 'ibd_files' in files_dict and files_dict['ibd_files'] is not None:
            files_dict['ibd_files'] = [Path(f) for f in files_dict['ibd_files']]
        if 'haproh_files' in files_dict and files_dict['haproh_files'] is not None:  # aDNA
            files_dict['haproh_files'] = [Path(f) for f in files_dict['haproh_files']]
        
        # Create nested config objects
        files_config = FilesConfig(**files_dict)
        algorithm_config = AlgorithmConfig(**config_dict.get("algorithm", {}))
        output_config = OutputConfig(**config_dict.get("output", {}))
        
        return cls(
            files=files_config,
            algorithm=algorithm_config, 
            output=output_config
        )

    @classmethod
    def from_cli_and_yaml(cls, cli_args: Dict[str, Any]) -> "PonderosaConfig":
        """Load from YAML file if provided, then override with CLI args."""
        
        # Start with YAML config if provided
        if "config" in cli_args and cli_args["config"]:
            config = cls.from_yaml(cli_args["config"])
            config_dict = config.to_dict()
        else:
            config_dict = {"files": {}, "algorithm": {}, "output": {}}
        
        # Override with CLI arguments
        # Map flat CLI args to nested structure
        file_args = [
            "ibd", "fam", "ages", "mapf", "populations", "training", "ibd_caller",
            # aDNA:
            "kin_file", "read_file"
        ]
        algorithm_args = [
            "min_segment_length", "min_total_ibd", "population", "genome_length",
            # aDNA:
            "pseudo_haploid", "min_overlap_snps", "trim_damage_bases",
            "min_pairs_for_normalization", "allow_no_segments"
        ]
        output_args = ["output", "min_probability", "verbose"]
        
        for arg, value in cli_args.items():
            if value is not None:
                if arg in file_args:
                    config_dict["files"][arg] = value
                elif arg in algorithm_args:
                    config_dict["algorithm"][arg] = value
                elif arg in output_args:
                    config_dict["output"][arg] = value
        
        return cls.from_dict(config_dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        # aDNA: added explicit, accurate keys; removed non-existent 'king'/'map' keys from earlier draft
        return {
            "files": {
                "ibd": str(self.files.ibd) if self.files.ibd else None,
                "ibd_files": [str(p) for p in (self.files.ibd_files or [])] if self.files.ibd_files else None,
                "fam": str(self.files.fam),
                "ibd_caller": self.files.ibd_caller,
                "ages": str(self.files.ages) if self.files.ages else None,
                "mapf": str(self.files.mapf) if self.files.mapf else None,
                "map_files": [str(p) for p in (self.files.map_files or [])] if self.files.map_files else None,
                "populations": str(self.files.populations) if self.files.populations else None,
                "training": str(self.files.training) if self.files.training else None,
                "rel_tree": str(self.files.rel_tree) if self.files.rel_tree else None,
                "priors": str(self.files.priors) if self.files.priors else None,
                # aDNA:
                "kin_file": str(self.files.kin_file) if self.files.kin_file else None,
                "read_file": str(self.files.read_file) if self.files.read_file else None,
                "haproh_file": str(self.files.haproh_file) if self.files.haproh_file else None,
                "haproh_files": [str(p) for p in (self.files.haproh_files or [])] if self.files.haproh_files else None,
            },
            "algorithm": {
                "min_segment_length": self.algorithm.min_segment_length,
                "min_total_ibd": self.algorithm.min_total_ibd,
                "max_gap": self.algorithm.max_gap,
                "use_phase_correction": self.algorithm.use_phase_correction,
                "mean_phase_error_distance": self.algorithm.mean_phase_error_distance,
                "degree_threshold": self.algorithm.degree_threshold,
                "haplotype_threshold": self.algorithm.haplotype_threshold,
                "population": self.algorithm.population,
                "genome_length": self.algorithm.genome_length,
                "validate_pair_order": self.algorithm.validate_pair_order,
                "parallel_processing": self.algorithm.parallel_processing,
                # aDNA:
                "pseudo_haploid": self.algorithm.pseudo_haploid,
                "min_overlap_snps": self.algorithm.min_overlap_snps,
                "trim_damage_bases": self.algorithm.trim_damage_bases,
                "min_pairs_for_normalization": self.algorithm.min_pairs_for_normalization,
                "allow_no_segments": self.algorithm.allow_no_segments,
            },
            "output": {
                "output": self.output.output,
                "min_probability": self.output.min_probability,
                "write_readable": self.output.write_readable,
                "write_pickle": self.output.write_pickle,
                "write_detailed": self.output.write_detailed,
                "write_training": self.output.write_training,
                "create_plots": self.output.create_plots,
                "plot_format": self.output.plot_format,
                "plot_dpi": self.output.plot_dpi,
                "log_level": self.output.log_level,
                "verbose": self.output.verbose
            }
        }
    
    def validate(self) -> None:
        """Validate the entire configuration."""
        self.files.validate()
        
        # Validate algorithm parameters
        if self.algorithm.min_segment_length <= 0:
            raise ValueError("min_segment_length must be positive")
        
        if self.algorithm.min_total_ibd < 0:
            raise ValueError("min_total_ibd must be non-negative")
        
        if not 0 <= self.output.min_probability <= 1:
            raise ValueError("min_probability must be between 0 and 1")

        # aDNA: soft checks; keep permissive
        if self.algorithm.min_overlap_snps < 0:
            raise ValueError("min_overlap_snps must be >= 0")
        if self.algorithm.trim_damage_bases is not None and self.algorithm.trim_damage_bases < 0:
            raise ValueError("trim_damage_bases must be >= 0")
        if self.algorithm.min_pairs_for_normalization < 0:
            raise ValueError("min_pairs_for_normalization must be >= 0")
        
        logger.info("Configuration validation passed")
    
    def save_yaml(self, output_path: Union[str, Path]) -> None:
        """Save configuration to YAML file."""
        output_path = Path(output_path)
        
        with open(output_path, 'w') as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False, indent=2)
        
        logger.info(f"Configuration saved to {output_path}")
