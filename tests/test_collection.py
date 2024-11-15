import unittest
from wipe.modules.collection import (
    compile_results_checkm2,
    compile_results_prodigal,
    generate_coords,
)


class CollectionTests(unittest.TestCase):
    def test_compile_results_checkm2(self):
        obs = compile_results_checkm2(
            "/projects/greengenes2/gg2_genomes/linearized_test/G"
        )
        obs.to_csv(
            "./tests/tmp/checkm2.tsv", index=False, header=True, sep="\t"
        )

    def test_compile_results_prodigal(self):
        obs = compile_results_prodigal(
            "/projects/greengenes2/gg2_genomes/linearized_test/G"
        )
        obs.to_csv(
            "./tests/tmp/proteins.tsv", index=False, header=True, sep="\t"
        )

    def test_generate_coords(self):
        generate_coords(
            "/projects/greengenes2/gg2_genomes/linearized_test/G/964/307",
            "./tests/tmp/coords.txt.xz",
        )


if __name__ == "__main__":
    pass
    # unittest.main()
