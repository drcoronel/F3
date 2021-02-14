def parser_test():
    from src.horizon_utils import parse_horizons_3D, sort_horizons
    import pathlib
    import matplotlib.pyplot as plt

    file = "data/Horizons_3D.dat"
    path = pathlib.Path(file)

    xyt, surfaces = parse_horizons_3D(path)

    horizon_dict = sort_horizons(xyt, surfaces)

    for k in horizon_dict.keys():
        temp = horizon_dict[k]
        plt.scatter(temp[0], temp[1], c=temp[2], cmap="viridis")
        plt.title(k)
        plt.show()
