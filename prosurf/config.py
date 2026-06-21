from dataclasses import dataclass, field
from pathlib import Path

@dataclass
class MetricConfig:
    rsasa_threshold: float = 0.20
    patch_radius_frac: float = 0.55  # radius = max(frac * R_g of charge cloud, 8 Å)
    his_weight: float = 0.0
    z_percentile: float = 90.0
    seed: int = 0
    adjacency_radius: float = 8.0

@dataclass
class PathsConfig:
    data_dir: Path = field(default_factory=lambda: Path("data"))
    cache_dir: Path = field(default_factory=lambda: Path("data/cache"))
