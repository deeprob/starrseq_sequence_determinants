import os
import utils as ut

DESC = "for peaks called in all libraries by starrpeaker, get the homer TF motifs that might bind to those peaks"


def main(meta_peak_file, homer_motif_file, genome, store_dir, method="homer", threads=64):

    # convert meta peak file to correct format
    formatted_peak_store_file = os.path.join(store_dir, f"{method}_in", "peaks.bed")
    os.makedirs(os.path.dirname(formatted_peak_store_file), exist_ok=True)
    ut.create_homer_compatible_bedfile(meta_peak_file, formatted_peak_store_file)

    # get database motif binding locations
    outdir = os.path.join(store_dir, f"{method}_out")
    os.makedirs(outdir, exist_ok=True)
    ut.pwm_scan_homer(formatted_peak_store_file, homer_motif_file, genome, outdir, threads)
    return


if __name__ == "__main__":
    meta_peak_file = "/data5/deepro/starrseq/main_library/4_quality_control_peaks/data/peak_cov/meta_peak_file.bed"
    homer_motif_file = "/data5/deepro/starrseq/main_library/6_sequence_determinants/data/homer/homer.motifs" 
    store_dir = "/data5/deepro/starrseq/main_library/6_sequence_determinants/data/meta_peaks"
    genome = "/data5/deepro/genomes/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    main(meta_peak_file, homer_motif_file, genome, store_dir)
