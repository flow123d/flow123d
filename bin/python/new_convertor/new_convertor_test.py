import unittest
import filecmp
from new_convertor import *


def files_cmp(ref,out):
    with open(ref, "r") as f:
        t=ruml.load(f, Loader=ruml.RoundTripLoader)
    with open(ref+".norm", "w") as f:
        ruml.dump(t, f, Dumper=ruml.RoundTripDumper)
    return filecmp.cmp(ref+".norm", out)

def remove_prefix(str, prefix):
    if str.startswith(prefix):
        return str[len(prefix):]
    return str

class TestActions(unittest.TestCase):

    def test_add_key(self):
        changes=Changes()
        changes.new_version("A")
        changes += PathSet(["/"]).add_key_to_map(key="flow123d_version", value="2.0.0")
        self.perform(changes)

    def perform(self, changes):
        in_file="test_input.yaml"
        with open(in_file, "r") as f:
            root = ruml.load(f, Loader=ruml.RoundTripLoader)

        changes.apply_changes([ (in_file, root) ],"A",None,reversed=False)

        test_name = self.id().split('.')[-1]
        test_name_base = remove_prefix(test_name, "test_")
        out_file = "test_"+test_name_base+".out.yaml"
        ref_file = "test_"+test_name_base+".ref.yaml"
        with open(out_file, "w") as f:
            ruml.dump(root, f, Dumper=ruml.RoundTripDumper)

        self.assertTrue(files_cmp(ref_file, out_file))

if __name__ == '__main__':
    unittest.main()