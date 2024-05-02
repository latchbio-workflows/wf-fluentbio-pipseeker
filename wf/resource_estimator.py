import os
from typing import Optional
from latch.types import LatchDir, LatchFile
import requests
from wf.configurations import GenomeType
from wf.configurations import get_mapping_reference


def get_s3_object_size(url):
    # Convert the url to a format usable by response.
    url = url.replace("s3://latch-public/", "https://latch-public.s3.amazonaws.com/")
    try:
        response = requests.head(url)
        if response.status_code == 200:
            size_in_bytes = response.headers.get('Content-Length')
            return int(size_in_bytes) if size_in_bytes is not None else None
        else:
            print(f"Failed to retrieve object: HTTP {response.status_code}")
    except requests.exceptions.RequestException as e:
        print(f"HTTP Request failed: {e}")
    return None


def get_downsample_factor(*, downsample_to, input_reads):
    # Requires both downsample_to and input_reads inputs.
    downsample_factor = None
    if downsample_to and input_reads:
        if input_reads != 0:
            downsample_factor = downsample_to / input_reads
    return downsample_factor


def get_fastqs_size_bytes(*, fastq_directory, downsample_factor):
    fastqs_size_bytes = 0
    for file in fastq_directory.iterdir():  # Using iterdir() to iterate over contents
        if isinstance(file, LatchFile) and file.path.endswith('.fastq.gz'):
            fastqs_size_bytes += file.size()
    fastqs_size_gb = fastqs_size_bytes / 1024 ** 3
    print(f'Input FASTQs size: {fastqs_size_gb:.2f} GB')
    if downsample_factor is not None:
        fastqs_size_gb *= downsample_factor

    return fastqs_size_bytes


def baseline_ram_estimator(*, fastqs_size_gb, num_threads):
    # Baseline RAM footprint estimation.
    # Coefficients are in units of GB; fastqs_size is in units of bytes
    # Convert bytes to GB before calculating; then reconvert to bytes.

    baseline_ram_bytes = 1024 ** 3 * (2.24 + (0.01 * fastqs_size_gb) +
                                      (0.0011 * num_threads * fastqs_size_gb))
    return baseline_ram_bytes


def barcoding_ram_estimator(*, baseline_ram_bytes, num_threads):
    # Barcoding memory estimation.
    # Coefficients are in units of GB; baseline_ram is in bytes; converting all but baseline_ram to bytes.
    # Since we can't easily measure the r1+r2 length here, we will default to a max expected value of 300bp.
    r1r2_length_sum = 300  # have to default to max of 300 since can't measure here.
    barcoding_ram = 1024 ** 3 * (
            (1.23 * num_threads) + (0.0166 * r1r2_length_sum) + (0.009 * num_threads * r1r2_length_sum))
    barcoding_ram_bytes = baseline_ram_bytes + barcoding_ram
    return barcoding_ram_bytes


def mapping_ref_size_estimator(*, genome_source, prebuilt_genome, custom_prebuilt_genome,
                               custom_prebuilt_genome_zipped):
    """
    Estimate the size of the mapping reference, depending on whether it is compressed or not.
    For simplicity, we take the size as 1.25x the size of the reference genome, which is equivalent to the
        filesize reduction using the maximum gzip compression level.
    """
    star_index_size_bytes = 0

    reference_p = get_mapping_reference(genome_source=genome_source, prebuilt_genome=prebuilt_genome,
                                        custom_prebuilt_genome=custom_prebuilt_genome,
                                        custom_prebuilt_genome_zipped=custom_prebuilt_genome_zipped,
                                        get_path_only=True)

    if genome_source == 'prebuilt_genome':
        # Prebuilt genomes are stored in latch s3 public bucket and must use special lookup for file size.
        #  Here, the reference_p is the s3_url rather than LatchFile instance.
        #  get_s3_object_size returns None if response fails to retrieve the size of the url.
        star_index_size_bytes = get_s3_object_size(reference_p)

        if star_index_size_bytes is None:
            #  response failed to obtain the size for the reference.
            #   Set default value (max size prebuilt ref is human + mouse, 24GB) and print warning.
            star_index_size_bytes = 24 * 1024 ** 3
        return star_index_size_bytes

    # Check whether the following file extensions are present in the reference path.
    #   We're working with a PosixPath obj.
    if reference_p:
        if reference_p.suffixes:
            # Input is a file. Check suffixes to confirm is compressed.
            if any([ext == reference_p.suffixes[0] for ext in ['.tar', '.tar.gz', '.zip']]):
                # Returns the unpacked reference size estimate.
                star_index_size_bytes = 1.25 * os.path.getsize(reference_p)
            else:
                star_index_size_bytes = os.path.getsize(reference_p)
        else:
            # Input is a directory. Iterate over contents to get size.
            try:
                star_index_size_bytes = sum(
                    [os.path.getsize(os.path.join(reference_p, p)) for p in os.listdir(reference_p)])
            except:
                ValueError(f"Could not calculate the size of the reference genome directory at {reference_p}. \n"
                           f"Folder contents: {os.listdir(reference_p)}")
    else:
        print('No reference path was found. Continuing without reference genome size estimation.')

    print(f'STAR index size: {star_index_size_bytes / 1024**3:.3f} GB  for {reference_p}')
    return star_index_size_bytes


def star_ram_estimator(*, star_index_size_bytes, num_threads, baseline_ram_bytes, sorted_bam, safety_margin=1.1):
    """
    Estimate STAR consumption with a 10% safety margin OR with 2x expected STAR RAM for when processing T100-1000s.

    Coefficients are in GB, so convert to bytes where needed;
    star_index_size is already in units of bytes.

    Params:
        safety_margin: Adds a 10% safety margin for case of STAR without counting.
    """
    # Baseline RAM requirement is added below.
    star_ram_bytes = (0.93 * star_index_size_bytes) + (0.55 + (0.23 * num_threads) * 1024 ** 3)

    if not sorted_bam:
        star_ram_bytes = safety_margin * (baseline_ram_bytes + star_ram_bytes)
    else:
        # Sorted BAM in STAR requires MI counting/correction --> T100 samples spike RAM by ~2x.
        star_ram_bytes = baseline_ram_bytes + 2 * star_ram_bytes
    return star_ram_bytes


def molinfo_ram_estimator(*, baseline_ram_bytes, fastqs_size_bytes, num_threads, exons_only=False):
    """
    Molecule info memory estimation (includes BAM-parsing estimation).
    Coefficients are in units of GB when not multiplied by bytes data.
    Convert all GB values to bytes data aside from baseline_ram and fastqs_size, which are already in bytes units.
    """
    fastqs_size_gb = fastqs_size_bytes / 1024 ** 3
    molinfo_ram_bytes = baseline_ram_bytes + (0.772 * fastqs_size_gb) + \
                        1024 ** 3 * ((0.54 * exons_only) + (0.525 * num_threads) +
                                     (0.71 * num_threads * exons_only))
    return molinfo_ram_bytes


def get_num_threads(fastq_directory: Optional[LatchDir] = None,
                    **kwargs) -> int:
    # Set number of cores based on a balance between cost savings and speed, based on the size of input fastqs.

    fastqs_size_bytes = get_fastqs_size_bytes(fastq_directory=fastq_directory, downsample_factor=None)
    fastqs_size_gb = fastqs_size_bytes / 1024 ** 3

    if fastqs_size_gb < 1:
        num_threads = 4
    elif fastqs_size_gb < 4:
        num_threads = 8
    elif fastqs_size_gb < 8:
        num_threads = 12
    elif fastqs_size_gb < 16:
        num_threads = 24
    elif fastqs_size_gb < 32:
        num_threads = 32
    else:
        num_threads = 64
    return num_threads


def get_disk_requirement_gb(*, fastq_directory: Optional[LatchDir] = None,
                            snt_fastq: Optional[LatchDir] = None,
                            hto_fastq: Optional[LatchDir] = None,
                            sorted_bam: bool = False,
                            downsample_to: Optional[int] = None,
                            input_reads: Optional[int] = None,
                            genome_source: str,
                            prebuilt_genome: GenomeType,
                            custom_prebuilt_genome: Optional[LatchDir],
                            custom_prebuilt_genome_zipped: Optional[LatchFile],
                            safety_margin=1.5,  # include 50% overage since have unexpected overhead.
                            **kwargs
                            ) -> int:
    # Set defaults.
    fastqs_size_bytes = 0
    snt_fastq_size_bytes = 0
    hto_fastq_size_bytes = 0

    if fastq_directory:
        downsample_factor = get_downsample_factor(downsample_to=downsample_to, input_reads=input_reads)
        fastqs_size_bytes = get_fastqs_size_bytes(fastq_directory=fastq_directory, downsample_factor=downsample_factor)

    if snt_fastq:
        snt_fastq_size_bytes = get_fastqs_size_bytes(fastq_directory=snt_fastq, downsample_factor=None)
    if hto_fastq:
        hto_fastq_size_bytes = get_fastqs_size_bytes(fastq_directory=hto_fastq, downsample_factor=None)

    if sorted_bam:
        # Previously 10x, but this option will output 2 BAM files (add 1x more space).
        #   Also add +1x for storing the input fastqs.
        fastqs_size_bytes = fastqs_size_bytes * 12 + snt_fastq_size_bytes + hto_fastq_size_bytes
    else:
        # 2.5x + 1x for housing the input fastqs
        fastqs_size_bytes = fastqs_size_bytes * 3.5 + snt_fastq_size_bytes + hto_fastq_size_bytes

    # Add on the STAR index size, since need to store and unpack in some cases.
    star_index_size_bytes = mapping_ref_size_estimator(genome_source=genome_source,
                                                        prebuilt_genome=prebuilt_genome,
                                                        custom_prebuilt_genome=custom_prebuilt_genome,
                                                        custom_prebuilt_genome_zipped=custom_prebuilt_genome_zipped)
    if prebuilt_genome is not None or custom_prebuilt_genome_zipped is not None:
        # If need to unpack the index, add on 2.25x the size of the index to cover the index + unpacked copy.
        star_index_size_bytes = star_index_size_bytes * 2.25

    # Final calculation
    required_space_gb = (fastqs_size_bytes + star_index_size_bytes + 1 ) * safety_margin / 1024 ** 3

    # In event have < 1GB, will be rounded to 0. Setting minimum default disk as 2GB.
    return int(max(required_space_gb, 2))


def get_memory_requirement_gb(*,
                              fastq_directory: Optional[LatchDir] = None,
                              genome_source: str,
                              prebuilt_genome: GenomeType,
                              custom_prebuilt_genome: Optional[LatchDir],
                              custom_prebuilt_genome_zipped: Optional[LatchFile],
                              downsample_to: Optional[int] = None,
                              input_reads: Optional[int] = None,
                              sorted_bam: bool = False,
                              **kwargs) -> int:
    # Get fastq size.
    downsample_factor = get_downsample_factor(downsample_to=downsample_to, input_reads=input_reads)
    fastqs_size_bytes = get_fastqs_size_bytes(fastq_directory=fastq_directory, downsample_factor=downsample_factor)
    fastqs_size_gb = fastqs_size_bytes / 1024 ** 3

    # Num threads calc.
    num_threads = get_num_threads(fastq_directory=fastq_directory)

    # Baseline ram calc.
    baseline_ram_bytes = baseline_ram_estimator(fastqs_size_gb=fastqs_size_gb, num_threads=num_threads)
    barcoding_ram_bytes = barcoding_ram_estimator(baseline_ram_bytes=baseline_ram_bytes, num_threads=num_threads)
    star_index_size_bytes = mapping_ref_size_estimator(genome_source=genome_source, prebuilt_genome=prebuilt_genome,
                                                       custom_prebuilt_genome=custom_prebuilt_genome,
                                                       custom_prebuilt_genome_zipped=custom_prebuilt_genome_zipped)
    star_ram_bytes = star_ram_estimator(star_index_size_bytes=star_index_size_bytes, num_threads=num_threads,
                                        baseline_ram_bytes=baseline_ram_bytes, sorted_bam=sorted_bam)
    molinfo_ram_bytes = molinfo_ram_estimator(baseline_ram_bytes=baseline_ram_bytes,
                                              fastqs_size_bytes=fastqs_size_bytes,
                                              num_threads=num_threads)

    return int(max(barcoding_ram_bytes, star_ram_bytes, molinfo_ram_bytes) / 1024 ** 3)
