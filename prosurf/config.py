from dataclasses import dataclass, field
from pathlib import Path

@dataclass
class MetricConfig:
    rsasa_threshold: float = 0.20
    patch_radius: float = 10.0
    his_weight: float = 0.0
    z_percentile: float = 90.0
    seed: int = 0

@dataclass
class PathsConfig:
    data_dir: Path = field(default_factory=lambda: Path("data"))
    cache_dir: Path = field(default_factory=lambda: Path("data/cache"))
