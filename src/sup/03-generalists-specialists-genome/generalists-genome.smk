rule genome_metadata:
    output:
        os.path.join("cache", CACHE_META, "{GENUS}.xml")
    threads: 1
    conda:
        "conda/ncbi-datasets-cli.yaml"
    shell:
        """
        datasets \
            summary genome \
            taxon {wildcards.GENUS} \
            --reference > {output}
        """

# import pandas as pd

# # Import generalists genus
# generalists = pd.read_csv("raw/generalists/generalist.txt", header = 0, names = ["generalist"], sep = "\t")
# genus = list(generalists["generalist"])
# genus = genus[0:2]

# # Cache directories
# CACHE_META = "cache/11-genome-metadata"

# # Functions to generate file names
# def get_genome_meta_files(wildcards):
#     global FILES
#     FILES = expand(os.path.join("cache", CACHE_META, "{GENUS}.xml"), GENUS = genus)
#     return FILES

# rule genome_metadata:
#     output:
#         os.path.join("cache", CACHE_META, "{GENUS}.xml")
#     conda:
#         "src/conda/ncbi-datasets-cli.yaml"
#     threads: 1
#     shell:
#         """
#         datasets \
#             summary genome \
#             taxon {wildcards.GENUS} \
#             --reference > {output}
#         """