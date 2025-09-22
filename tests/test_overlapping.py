import os
import sys

from satellome.core_functions.io.tab_file import sc_iter_tab_file
from satellome.core_functions.models.trf_model import TRModel
from satellome.trf_embedings import get_cosine_distance

# Use environment variable or command line argument for test data path
# Set TEST_DATA_PATH environment variable to point to your test data directory
# e.g., export TEST_DATA_PATH=/path/to/test/data
input_file = os.environ.get('TEST_DATA_PATH',
                           sys.argv[1] if len(sys.argv) > 1
                           else "test_data/GCF_905171775.1_aRanTem1.1_genomic.1kb.trf")
# Check if input file exists
if not os.path.exists(input_file):
    print(f"Error: Input file not found: {input_file}")
    print("Please provide a valid TRF file path via:")
    print("  - Command line argument: python test_overlapping.py /path/to/file.trf")
    print("  - Environment variable: export TEST_DATA_PATH=/path/to/file.trf")
    sys.exit(1)

last_trf = None
for j, trf_obj in enumerate(sc_iter_tab_file(input_file, TRModel)):
    if not last_trf:
        last_trf = trf_obj
        continue
    #     print(last_trf.trf_l_ind, last_trf.trf_r_ind, trf_obj.trf_l_ind, trf_obj.trf_r_ind)
    if (
        trf_obj.trf_head == last_trf.trf_head
        and last_trf.trf_r_ind > trf_obj.trf_l_ind
        and last_trf.trf_r_ind < trf_obj.trf_r_ind
    ):
        #         print(trf_obj)
        #         print(last_trf)
        vector1 = trf_obj.get_vector()
        vector2 = last_trf.get_vector()
        dist = get_cosine_distance(vector1, vector2)
        print(
            last_trf.trf_l_ind,
            last_trf.trf_r_ind,
            trf_obj.trf_l_ind,
            trf_obj.trf_r_ind,
            dist,
        )
        print(last_trf.trf_consensus, trf_obj.trf_consensus)
    #         break
    last_trf = trf_obj
