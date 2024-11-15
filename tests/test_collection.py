import unittest
from wipe.modules.collection import (
    compile_results_checkm2,
    compile_results_prodigal,
    generate_coords,
)


class CollectionTests(unittest.TestCase):
    def test_compile_results_checkm2(self):
        pass

    # def test_compile_results_prodigal(self):
    #     obs = compile_results_prodigal(
    #         "/projects/greengenes2/gg2_genomes/linearized_test/G"
    #     )
    #     print(obs.head(5))

    # def test_generate_coords(self):
    #     generate_coords(
    #         "/projects/greengenes2/gg2_genomes/linearized_test/G/964/307",
    #         "./tests/tmp/coords.txt.xz",
    #     )
