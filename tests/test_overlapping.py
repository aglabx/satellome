from satellome.core_functions.io.tab_file import sc_iter_tab_file
from satellome.core_functions.models.trf_model import TRModel
from satellome.trf_embedings import get_cosine_distance

input_file = "/mnt/data/podgornaya/rana_temporaria/users/akomissarov/trf/GCF_905171775.1_aRanTem1.1_genomic.1kb.trf"
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
