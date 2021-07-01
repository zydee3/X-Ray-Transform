from x_ray_transform import Surface, Euclidean


def test_surface():
    surface = Surface()
    metric = Euclidean()
    surface.set_metric_as_euclidean(metric)
