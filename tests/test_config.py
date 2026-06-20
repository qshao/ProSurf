from prosurf.config import MetricConfig

def test_metric_config_defaults():
    cfg = MetricConfig()
    assert cfg.rsasa_threshold == 0.20
    assert cfg.patch_radius == 10.0
    assert cfg.his_weight == 0.0
    assert cfg.z_percentile == 90.0
    assert cfg.seed == 0

def test_metric_config_override():
    cfg = MetricConfig(his_weight=0.1)
    assert cfg.his_weight == 0.1
    assert cfg.rsasa_threshold == 0.20
