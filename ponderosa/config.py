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

    # aDNA: allow optional ROH inputs (e.g., hapROH output or precomputed ROH bed/tsv)
    roh: Optional[Path] = None                     # single ROH file
    roh_files: Optional[List[Path]] = None         # list of per-chrom or multi-file ROH
    roh_format: Optional[str] = None               # e.g., "haproh", "bed", "eigenstrat-roh"
    _roh_file_list: Optional[List[Path]] = field(default=None, init=False)

    # aDNA: optional KING-like kinship table to seed founders/sim (used elsewhere, keep here for consistency)
    king: Optional[Path] = None

    def _validate_file_lists(self, single_file: Path, file_list: List[Path], file_type: str) -> List[Path]:

        single_file_arg = lambda x: {"map": "mapf", "ibd": "ibd", "roh": "roh"}[x]  # aDNA: added roh
        multiple_file_arg = lambda x: f"{x}_files"

        if single_file and file_list:
            raise ValueError(f"Cannot specify both '{single_file_arg(file_type)}' and '{multiple_file_arg(file_type)}'. Choose one.")

        if single_file:
            if not single_file.exists():
                raise FileNotFoundError(f"{single_file_arg(file_type)} file not found: {single_file}")
            
            # Preserve original behavior for .map pattern completion
            if file_type == "map" and "chr1" in single_file.name: # Populate with all other chromosomes
                out_list = [single_file.with_name(single_file.name.replace("chr1", f"chr{i}")) for i in range(1, 23)]
            else:
                out_list = [single_file]

        elif file_list:
            if not file_list:
                raise ValueError(f"{multiple_file_arg(file_type)} list cannot be empty")
            out_list = []
            for i, file_n in enumerate(file_list):
                if not Path(file_n).exists():
                    raise FileNotFoundError(f"{multiple_file_arg(file_type)} file not found: {file_n} (index {i})")
                out_list.append(Path(file_n))

        else:
            # NOTE: original code had a ValueError but didn’t raise it — fix to actually raise.
            raise ValueError(f"{file_type} file(s) have not been provided. "
                             f"Specify either '{single_file_arg(file_type)}' or '{multiple_file_arg(file_type)}'.")

        return out_list
    
    def validate(self) -> None:
        """Validate that required files exist."""
        # Check required files
        required_files = [self.fam]
        for file_path in required_files:
            if not file_path.exists():
                raise FileNotFoundError(f"Required file not found: {file_path}")
            
        # map + ibd lists (existing behavior)
        if self.ibd or self.ibd_files:
            self._ibd_file_list = self._validate_file_lists(self.ibd, self.ibd_files, "ibd")
        else:
            # Keep optional (some runs may be training-only or simulation)
            self._ibd_file_list = []

        if self.mapf or self.map_files:
            self._map_file_list = self._validate_file_lists(self.mapf, self.map_files, "map")
        else:
            self._map_file_list = []

        # aDNA: validate optional roh lists if provided
        if self.roh or self.roh_files:
            self._roh_file_list = self._validate_file_lists(self.roh, self.roh_files, "roh")
        else:
            self._roh_file_list = []

        # Check optional files if provided
        optional_files = [self.ages, self.mapf, self.populations, self.training, self.priors, self.king, self.roh]
        for file_path in optional_files:
            if file_path and not Path(file_path).exists():
                logger.warning(f"Optional file not found: {file_path}")

        # aDNA: light validation of roh_format when ROH files exist
        if self._roh_file_list and not self.roh_format:
            # default to hapROH when ROH files are provided but no format specified
            self.roh_format = "haproh"
        if self.roh_format:
            allowed = {"haproh", "bed", "eigenstrat-roh"}
            if self.roh_format.lower() not in allowed:
                raise ValueError(f"Invalid roh_format='{self.roh_format}'. Allowed: {sorted(allowed)}")

    @property
    def map_file_list(self) -> List[Path]:
        """Return unified list of map files. Must call validate() first."""
        if self._map_file_list is None:
            raise RuntimeError("Must call validate() before accessing map_file_list")
        return self._map_file_list
    
    @property
    def ibd_file_list(self) -> List[Path]:
        """Return unified list of IBD files. Must call validate() first."""
        if self._ibd_file_list is None:
            raise RuntimeError("Must call validate() before accessing ibd_file_list")
        return self._ibd_file_list

    @property
    def roh_file_list(self) -> List[Path]:
        """aDNA: unified list of ROH files (hapROH/bed/etc). Must call validate() first."""
        if self._roh_file_list is None:
            raise RuntimeError("Must call validate() before accessing roh_file_list")
        return self._roh_file_list


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

    # aDNA: optional knobs to sensibly integrate ROH and ancient data quality (kept OFF by default)
    ancient_mode: bool = False                        # enable aDNA defaults elsewhere (no behavior change here)
    use_roh_priors: bool = False                      # when True, allow ROH-informed priors downstream
    roh_min_cm: float = 4.0                           # minimum ROH length (cM) to consider in summaries
    roh_merge_gap_cm: float = 0.5                     # merge gap for ROH stitching in summaries
    use_damage_aware: bool = False                    # if pipeline supports PMD-aware genotype handling
    min_coverage: float = 0.0                         # optional per-sample min coverage recorded by loaders
    max_contamination: float = 1.0                    # optional per-sample max contamination tolerated (record only)


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

    # aDNA: (optional) write ROH summaries
    write_roh_summary: bool = False                  # when True, downstream can emit {prefix}_roh_summary.txt


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
        # aDNA: include new optional paths ('king', 'roh')
        single_path_fields = {'ibd', 'fam', 'ages', 'mapf', 'populations', 'training', 'rel_tree', 'king', 'roh'}
        for key, value in files_dict.items():
            if value is not None and key in single_path_fields:
                files_dict[key] = Path(value)
        
        # Handle list fields (convert each string to Path)
        if 'map_files' in files_dict and files_dict['map_files'] is not None:
            files_dict['map_files'] = [Path(f) for f in files_dict['map_files']]
        if 'ibd_files' in files_dict and files_dict['ibd_files'] is not None:
            files_dict['ibd_files'] = [Path(f) for f in files_dict['ibd_files']]
        if 'roh_files' in files_dict and files_dict['roh_files'] is not None:     # aDNA
            files_dict['roh_files'] = [Path(f) for f in files_dict['roh_files']]
        
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
            "ibd", "fam", "ages", "mapf", "populations", "training",
            "ibd_caller", "rel_tree", "king", "roh"  # aDNA
        ]
        algorithm_args = [
            "min_segment_length", "min_total_ibd", "population", "genome_length",
            # aDNA knobs (optional)
            "ancient_mode", "use_roh_priors", "roh_min_cm", "roh_merge_gap_cm",
            "use_damage_aware", "min_coverage", "max_contamination"
        ]
        output_args = ["output", "min_probability", "verbose", "write_roh_summary"]  # aDNA
        
        for arg, value in cli_args.items():
            if value is not None:
                if arg in file_args:
                    config_dict["files"][arg] = value
                elif arg in algorithm_args:
                    config_dict["algorithm"][arg] = value
                elif arg in output_args:
                    config_dict["output"][arg] = value
                elif arg in ("map_files", "ibd_files", "roh_files"):    # aDNA lists
                    config_dict.setdefault("files", {})[arg] = value
        
        return cls.from_dict(config_dict)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        # NOTE: original version referenced self.files.king and self.files.map
        # while FilesConfig had no 'king' and used 'mapf'. We added 'king' above and
        # keep the outward YAML key as 'map' for compatibility.
        return {
            "files": {
                "ibd": str(self.files.ibd) if self.files.ibd else None,
                "fam": str(self.files.fam),
                "king": str(self.files.king) if self.files.king else None,            # aDNA: present but optional
                "ages": str(self.files.ages) if self.files.ages else None,
                "map": str(self.files.mapf) if self.files.mapf else None,             # keep outward name 'map'
                "map_files": [str(p) for p in (self.files.map_files or [])] or None,
                "ibd_files": [str(p) for p in (self.files.ibd_files or [])] or None,
                "populations": str(self.files.populations) if self.files.populations else None,
                "training": str(self.files.training) if self.files.training else None,
                "rel_tree": str(self.files.rel_tree) if self.files.rel_tree else None,
                # aDNA outputs
                "roh": str(self.files.roh) if self.files.roh else None,
                "roh_files": [str(p) for p in (self.files.roh_files or [])] or None,
                "roh_format": self.files.roh_format
            },
            "algorithm": {
                "min_segment_length": self.algorithm.min_segment_length,
                "min_total_ibd": self.algorithm.min_total_ibd,
                "max_gap": self.algorithm.max_gap,
                "use_phase_correction": self.algorithm.use_phase_correction,
                "population": self.algorithm.population,
                "genome_length": self.algorithm.genome_length,
                "validate_pair_order": self.algorithm.validate_pair_order,
                "parallel_processing": self.algorithm.parallel_processing,
                # aDNA knobs
                "ancient_mode": self.algorithm.ancient_mode,
                "use_roh_priors": self.algorithm.use_roh_priors,
                "roh_min_cm": self.algorithm.roh_min_cm,
                "roh_merge_gap_cm": self.algorithm.roh_merge_gap_cm,
                "use_damage_aware": self.algorithm.use_damage_aware,
                "min_coverage": self.algorithm.min_coverage,
                "max_contamination": self.algorithm.max_contamination,
            },
            "output": {
                "output": self.output.output,
                "min_probability": self.output.min_probability,
                "write_readable": self.output.write_readable,
                "verbose": self.output.verbose,
                "write_training": self.output.write_training,
                "write_pickle": self.output.write_pickle,
                "write_detailed": self.output.write_detailed,
                # aDNA
                "write_roh_summary": self.output.write_roh_summary,
                "create_plots": self.output.create_plots,
                "plot_format": self.output.plot_format,
                "plot_dpi": self.output.plot_dpi,
                "log_level": self.output.log_level
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
        
        # aDNA: sanity checks (do not enforce behavior yet)
        if self.algorithm.roh_min_cm < 0:
            raise ValueError("roh_min_cm must be non-negative")
        if self.algorithm.roh_merge_gap_cm < 0:
            raise ValueError("roh_merge_gap_cm must be non-negative")
        if not (0.0 <= self.algorithm.max_contamination <= 1.0):
            raise ValueError("max_contamination must be within [0,1]")

        logger.info("Configuration validation passed")
    
    def save_yaml(self, output_path: Union[str, Path]) -> None:
        """Save configuration to YAML file."""
        output_path = Path(output_path)
        
        with open(output_path, 'w') as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False, indent=2)
        
        logger.info(f"Configuration saved to {output_path}")

