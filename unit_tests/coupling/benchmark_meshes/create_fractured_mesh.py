import os, sys
import yaml
import numpy as np
from typing import Tuple, List
from endorse import common

from bgem.gmsh import gmsh, options, gmsh_io
from bgem.stochastic import fracture

script_dir = os.path.dirname(os.path.realpath(__file__))


def box_with_sides(factory, dimensions):
    """
    Make a box and dictionary of its sides named: 'side_[xyz][01]'
    :return: box, sides_dict
    """
    box = factory.box(dimensions).set_region("box")
    side_z = factory.rectangle([dimensions[0], dimensions[1]])
    side_y = factory.rectangle([dimensions[0], dimensions[2]])
    side_x = factory.rectangle([dimensions[2], dimensions[1]])
    sides = dict(
        side_z0=side_z.copy().translate([0, 0, -dimensions[2] / 2]),
        side_z1=side_z.copy().translate([0, 0, +dimensions[2] / 2]),
        side_y0=side_y.copy().translate([0, 0, -dimensions[1] / 2]).rotate([-1, 0, 0], np.pi / 2),
        side_y1=side_y.copy().translate([0, 0, +dimensions[1] / 2]).rotate([-1, 0, 0], np.pi / 2),
        side_x0=side_x.copy().translate([0, 0, -dimensions[0] / 2]).rotate([0, 1, 0], np.pi / 2),
        side_x1=side_x.copy().translate([0, 0, +dimensions[0] / 2]).rotate([0, 1, 0], np.pi / 2)
    )
    for name, side in sides.items():
        side.modify_regions(name)
    return box, sides

def create_fractures_rectangles(gmsh_geom, fractures, shift, base_shape: 'ObjectSet'):
    # From given fracture date list 'fractures'.
    # transform the base_shape to fracture objects
    # fragment fractures by their intersections
    # return dict: fracture.region -> GMSHobject with corresponding fracture fragments
    if len(fractures) == 0:
        return []


    shapes = []
    for i, fr in enumerate(fractures):
        shape = base_shape.copy()
        print("fr: ", i, "tag: ", shape.dim_tags)
        shape = shape.scale([fr.rx, fr.ry, 1]) \
            .rotate(axis=[0,0,1], angle=fr.shape_angle) \
            .rotate(axis=fr.rotation_axis, angle=fr.rotation_angle) \
            .translate(fr.center + shift).set_region(fr.region)

        shapes.append(shape)

    fracture_fragments = gmsh_geom.fragment(*shapes)
    return fracture_fragments

def fr_dict_repr(fr):
    return dict(r=float(fr.r), normal=fr.normal.tolist(), center=fr.center.tolist(),
                aspect=float(fr.aspect), shape_angle=float(fr.shape_angle), region=fr.region.name)


def generate_fractures(pop:fracture.Population, range: Tuple[float, float], fr_limit, box,  seed) -> List[fracture.Fracture]:
    """
    Generate set of stochastic fractures.
    """
    np.random.seed(seed)
    max_fr_size = np.max(box)
    r_min, r_max = range
    if r_max is None:
        r_max = max_fr_size
    if r_min is None:
        # smallest size range
        n_frac_lim = fr_limit
    else:
        # prescribed fracture range
        n_frac_lim = None
    pop.domain = [b if d > 0 else 0.0 for d, b in zip(pop.domain, box)]
    pop.set_sample_range([r_min, r_max], sample_size=n_frac_lim)
    # logging.info(f"fr set range: {[r_min, r_max]}, fr_lim: {n_frac_lim}, mean population size: {pop.mean_size()}")


    pos_gen = fracture.UniformBoxPosition(pop.domain)
    fractures = pop.sample(pos_distr=pos_gen, keep_nonempty=True)
    for i, fr in enumerate(fractures):
        reg = gmsh.Region.get(f"fr_{i}")
        fr.region = reg

    # fracture.fr_intersect(fractures)
    return fractures


def fracture_set(cfg, fr_population:fracture.Population, seed):
    main_box_dimensions = cfg.geometry.box_dimensions

    # Fixed large fractures
    fix_seed = cfg.fractures.fixed_seed
    large_min_r = cfg.fractures.large_min_r
    large_box_dimensions = cfg.fractures.large_box
    fr_limit = cfg.fractures.n_frac_limit
    # logging.info(f"Large fracture seed: {fix_seed}")
    max_large_size = max([fam.size.diam_range[1] for fam in fr_population.families])
    fractures = generate_fractures(fr_population, (large_min_r, max_large_size), fr_limit, large_box_dimensions, fix_seed)

    large_fr_dict=dict(seed=fix_seed, fr_set=[fr_dict_repr(fr) for fr in fractures])
    with open(f"large_Fr_set.yaml", "w") as f:
        yaml.dump(large_fr_dict, f, sort_keys=False)
    n_large = len(fractures)
    #if n_large == 0:
    #    raise ValueError()
    # random small scale fractures
    small_fr = generate_fractures(fr_population, (None, large_min_r), fr_limit, main_box_dimensions, seed)
    fractures.extend(small_fr)
    # logging.info(f"Generated fractures: {n_large} large, {len(small_fr)} small.")
    return fractures, n_large


def make_geometry(factory, cfg_geom:'dotdict', cfg_mesh:'dotdict', fractures):
    box, sides = box_with_sides(factory, cfg_geom.box_dimensions)

    fractures = create_fractures_rectangles(factory, fractures, [0,0,0], factory.rectangle())
    fractures_group = factory.group(*fractures).intersect(box)

    box_fr, fractures_fr = factory.fragment(box, fractures_group)
    fractures_fr.mesh_step(cfg_mesh.fracture_mesh_step)  # .set_region("fractures")

    b_box_fr = box_fr.get_boundary().split_by_dimension()[2]
    b_fractures_fr = fractures_fr.get_boundary().split_by_dimension()[1]

    # select outer boundary
    boundary_mesh_step = cfg_mesh.boundary_mesh_step
    bbb_box = box_fr.get_boundary().copy()
    b_box = b_box_fr.select_by_intersect(bbb_box).set_region(".box_outer").mesh_step(
        boundary_mesh_step)
    b_fractures = b_fractures_fr.select_by_intersect(bbb_box).set_region(".fr_outer").mesh_step(
        boundary_mesh_step)

    boundary = factory.group(b_box, b_fractures)
    bulk_geom = factory.group(box_fr, fractures_fr, boundary)

    # without fractures
    # b_box = box.get_boundary()
    # boundary = factory.group(b_box)
    # bulk_geom = factory.group(box, boundary)

    return bulk_geom


def meshing(factory, objects, mesh_filename):
    """
    Common EDZ and transport domain meshing setup.
    """
    # factory.write_brep()
    #factory.mesh_options.CharacteristicLengthMin = cfg.get("min_mesh_step", cfg.boreholes_mesh_step)
    #factory.mesh_options.CharacteristicLengthMax = cfg.boundary_mesh_step
    factory.mesh_options.MinimumCirclePoints = 6
    factory.mesh_options.MinimumCurvePoints = 3
    #factory.mesh_options.Algorithm = options.Algorithm3d.MMG3D

    # mesh.Algorithm = options.Algorithm2d.MeshAdapt # produce some degenerated 2d elements on fracture boundaries ??
    # mesh.Algorithm = options.Algorithm2d.Delaunay
    # mesh.Algorithm = options.Algorithm2d.FrontalDelaunay

    factory.mesh_options.Algorithm = options.Algorithm3d.Delaunay
    #mesh.ToleranceInitialDelaunay = 0.01
    # mesh.ToleranceEdgeLength = fracture_mesh_step / 5
    #mesh.CharacteristicLengthFromPoints = True
    #factory.mesh_options.CharacteristicLengthFromCurvature = False
    #factory.mesh_options.CharacteristicLengthExtendFromBoundary = 2  # co se stane if 1
    #mesh.CharacteristicLengthMin = min_el_size
    #mesh.CharacteristicLengthMax = max_el_size

    #factory.keep_only(*objects)
    #factory.remove_duplicate_entities()
    factory.make_mesh(objects, dim=3)
    factory.write_mesh(filename=mesh_filename, format=gmsh.MeshFormat.msh2)


def make_gmsh(fractures:List['Fracture'], cfg:'dotdict'):
    """
    :param cfg_geom: repository mesh configuration cfg.repository_mesh
    :param fractures:  generated fractures
    :param mesh_file:
    :return:
    """
    mesh_filename =  os.path.join(cfg.output_dir, cfg.mesh_name + ".msh2")
    factory = gmsh.GeometryOCC(cfg.mesh_name, verbose=True)
    factory.get_logger().start()
    gopt = options.Geometry()
    gopt.Tolerance = 0.0001
    gopt.ToleranceBoolean = 0.001

    bulk = make_geometry(factory, cfg.geometry, cfg.mesh, fractures)

    meshing(factory, [bulk], mesh_filename)
    # factory.show()
    del factory
    return common.File(mesh_filename)


def make_mesh(workdir, output_dir, cfg_file):
    conf_file = os.path.join(workdir, cfg_file)
    cfg = common.config.load_config(conf_file)
    cfg.output_dir = output_dir

    fr_population = fracture.Population.initialize_3d(cfg.fractures.population, cfg.geometry.box_dimensions)
    fractures, n_large = fracture_set(cfg, fr_population, cfg.seed)

    mesh_file = make_gmsh(fractures, cfg)

    # the number of elements written by factory logger does not correspond to actual count
    reader = gmsh_io.GmshIO(mesh_file.path)
    print("N Elements: ", len(reader.elements))

    print("Mesh file: ", mesh_file)


if __name__ == '__main__':
    output_dir = None
    len_argv = len(sys.argv)
    assert len_argv > 1, "Specify input yaml file and output dir!"
    if len_argv == 2:
        output_dir = os.path.abspath(sys.argv[1])

    make_mesh(script_dir, output_dir, "./cube_fr_mesh_7f_3k.yaml")
    make_mesh(script_dir, output_dir, "./cube_fr_mesh_70f_30k.yaml")
    make_mesh(script_dir, output_dir, "./cube_fr_mesh_133f_300k.yaml")

