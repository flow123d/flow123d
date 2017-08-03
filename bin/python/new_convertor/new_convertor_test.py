import unittest
import filecmp
from shutil import copyfile
from new_convertor import *
import sys

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
logging.debug("First debug message.")
def files_cmp(ref,out):
    with open(ref, "r") as f:
        t=ruml.load(f, Loader=ruml.RoundTripLoader)
    with open(ref, "w") as f:
        ruml.dump(t, f, Dumper=ruml.RoundTripDumper)
    return filecmp.cmp(ref, out)

def remove_prefix(str, prefix):
    if str.startswith(prefix):
        return str[len(prefix):]
    return str


class TestActions(unittest.TestCase):

    def test_add_key(self):
        changes=Changes()
        changes.new_version("A")
        changes.add_key_to_map(PathSet(["/"]), key="flow123d_version", value="2.0.0")
        self.perform(changes)


    def test_rename_key(self):
        changes=Changes()
        changes.new_version("A")
        # Change degree keys in PadeApproximant
        path_set = PathSet(["/problem/secondary_equation/**/ode_solver!PadeApproximant/"])
        changes.rename_key(path_set, old_key="denominator_degree", new_key="pade_denominator_degree")
        changes.rename_key(path_set, old_key="nominator_degree", new_key="pade_nominator_degree")
        self.perform(changes)

    def test_rename_tag(self):
        changes=Changes()
        changes.new_version("A")

        # Rename equations and couplings
        path = PathSet(["/problem/secondary_equation/"])
        changes.rename_tag(path, old_tag="TransportOperatorSplitting", new_tag="Coupling_OperatorSplitting")
        self.perform(changes)

    """
    def test_replace_value(self):
        changes=Changes()
        changes.new_version("A")


        # Change sign of boundary fluxes
        path_set = PathSet(["/problem/(primary_equation|secondary_equation)/input_fields/*/bc_flux!FieldFormula/value/"])

        changes.replace_value(path_set,
            re_forward=('^(.*)$', '-(\\1)'),
            re_backward=('^(.*)$', '-(\\1)'))

        # Test manual conversion, alternative paths
        path_set = PathSet(["/problem/(secondary_equation|primary_equation)/input_fields/*/bc_type/"])

        changes.replace_value(path_set,
            re_forward=("^(robin|neumann)$", "total_flux"),
            re_backward=(None, "Select either 'robin' or 'neumann' according to the value of 'bc_flux', 'bc_pressure', 'bc_sigma'."))
        self.perform(changes)


    def test_scale_scalar(self):
        changes=Changes()
        changes.new_version("A")

        # Change sign of boundary fluxes
        path_set = PathSet([
            "/problem/(primary|secondary)_equation/input_fields/*/bc_flux/",
            "/problem/(primary|secondary)_equation/input_fields/*/bc_flux/#/",
            "/problem/(primary|secondary)_equation/input_fields/*/bc_flux!FieldConstant/value/",
        ])
        changes.scale_scalar(path_set, multiplicator=-1)
        self.perform(changes)
    """

####################################

    def  make_test_file(self, ext):
        fname = "tests/" + self.test_name_base + ext
        if not os.path.isfile(fname) and hasattr(self, "in_file"):
            copyfile(self.in_file, fname)
        return fname


    def perform(self, changes):
        test_name = self.id().split('.')[-1]
        self.test_name_base = remove_prefix(test_name, "test_")
        in_file = self.in_file = self.make_test_file(".in.yaml")
        out_file = self.make_test_file(".out.yaml")
        ref_file = self.make_test_file(".ref.yaml")
        rev_file = self.make_test_file(".rev.yaml")
        rrf_file = self.make_test_file(".rrf.yaml")

        with open(in_file, "r") as f:
            root = yml.load(f)
        changes.apply_changes([ (in_file, root) ],"A",None,reversed=False)
        with open(out_file, "w") as f:
            yml.dump(root, f)
        self.assertTrue(files_cmp(ref_file, out_file))

        with open(ref_file, "r") as f:
            root = yml.load(f)
        changes.apply_changes([ (ref_file, root) ],"A",None,reversed=True)
        with open(rev_file, "w") as f:
            yml.dump(root, f)
        self.assertTrue(files_cmp(rrf_file, rev_file))


if __name__ == '__main__':
    unittest.main()