# from app import utils
# from constants import *
# import app

# if __name__ == "__main__":
#     # Replace this with relevant test cases for your method
#     result = utils.generateBindingPredictions("6fa3cf66-72a7-4e9d-9165-a4a46cfd38c0", "HLA-A0201", Class_One_Predictors.MHCflurry)
#     print(result)

# from ProtPeptigram.DataProcessor import PeptideDataProcessor
# from ProtPeptigram.viz import ImmunoViz

# # Initialize data processor
# processor = PeptideDataProcessor()

# # Load data
# processor.load_peaks_data("data/peptides.csv")
# processor.load_protein_sequences("data/proteome.fasta")

# # Process data
# formatted_data = processor.filter_and_format_data(
#     filter_contaminants=True,
#     intensity_threshold=1000,
#     min_samples=2
# )

# # Create visualizations
# viz = ImmunoViz(formatted_data)
# fig, _ = viz.plot_peptigram(
#     protein_ids=["P20152", "P32261"],
#     group_by="Sample",
#     color_by="protein",
#     title="HLA Peptide Visualization"
# )

# # Save visualization
# fig.savefig("protein_visualization.png", dpi=300, bbox_inches="tight")

from unittest import result
from app import utils
from constants import *
import os
project_root = os.path.dirname(os.path.realpath(os.path.join(__file__, "..")))


if __name__ == "__main__":
    # Example 1 — no restricted allele list
    utils.run_clust_search(
        input_file="/mnt/c/Users/pmun0004/repos/Immunolyser/app/static/images/9f5c6456-f2fc-4ce2-889d-88b89ab6a1fd/GC_P7089/gibbscluster/P7089",
        ref_file=f"{project_root}/Immunolyser/app/tools/HLA-PepClust/data/ref_data",
        output_dir="/mnt/c/Users/pmun0004/repos/Immunolyser/app/static/images/9f5c6456-f2fc-4ce2-889d-88b89ab6a1fd/GC_P7089/hla_clust_output/P7089",
        species="human",
        # use_mhc_tp_full_DB="yes",   # or None — not "no"
        allele_file=os.path.join(project_root, 'app', 'static', 'mhc-tp-default-search-alleles - top 20.csv'),
        logger=None
    )

    # # Example 2 — use a restricted allele CSV (will append -hla ...)
    # utils.run_clust_search(
    #     input_file="/mnt/c/Users/pmun0004/repos/Immunolyser/app/static/images/9f5c6456-f2fc-4ce2-889d-88b89ab6a1fd/GC_P7089/gibbscluster/P7089",
    #     ref_file=f"{project_root}/app/tools/HLA-PepClust/data/ref_data",
    #     output_dir="/mnt/c/Users/pmun0004/repos/Immunolyser/app/static/images/9f5c6456-f2fc-4ce2-889d-88b89ab6a1fd/GC_P7089/hla_clust_output/P7089",
    #     species="human",
    #     use_mhc_tp_full_DB="no",                       # must be "no" (string) to enable restriction
    #     allele_file="/path/to/alleles.csv",            # CSV must contain 'Allele name standardised' header
    #     logger=None
    # )
