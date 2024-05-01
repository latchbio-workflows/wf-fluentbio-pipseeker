from enum import Enum
import subprocess
from pathlib import Path
from latch.types import LatchFile


class GenomeType(Enum):
    human = "Human"
    mouse = "Mouse"
    human_mouse = "Human and Mouse"
    drosophilia = "Drosophilia"
    zebrafish = "Zebrafish"
    arabidopsis_thaliana = "Arabidopsis thaliana"


class Chemistry(Enum):
    v3 = "v3"
    v4 = "v4"
    v5 = "v5"


class Verbosity(Enum):
    zero = "0"
    one = "1"
    two = "2"


def get_mapping_reference(*, genome_source, prebuilt_genome, custom_prebuilt_genome,
                          custom_prebuilt_genome_zipped, get_path_only=False):
    """
    Download and unpack prebuilt mapping references (using Fluent-built or custom).

    Args:
        get_path_only: Return only the compressed file path for filesize estimation (dynamic resource allocation).

    Returns:
        reference_p: Path to the prebuilt mapping reference.
    """
    reference_zipped_p = None

    print("\nPreparing reference genome")

    if genome_source == "prebuilt_genome":
        if prebuilt_genome == GenomeType.human:
            reference_zipped_p = LatchFile(
                "s3://latch-public/test-data/18440/pipseeker-gex-reference-GRCh38-2022.04.tar.gz"
            ).local_path
            reference_p = Path("/root/pipseeker-gex-reference-GRCh38-2022.04")

        elif prebuilt_genome == GenomeType.mouse:
            reference_zipped_p = LatchFile(
                "s3://latch-public/test-data/18440/pipseeker-gex-reference-GRCm39-2022.04.tar.gz"
            ).local_path
            reference_p = Path("/root/pipseeker-gex-reference-GRCm39-2022.04")

        elif prebuilt_genome == GenomeType.human_mouse:
            reference_zipped_p = LatchFile(
                "s3://latch-public/test-data/18440/pipseeker-gex-reference-GRCh38-and-GRCm39-2022.04.tar.gz"
            ).local_path
            reference_p = Path(
                "/root/pipseeker-gex-reference-GRCh38-and-GRCm39-2022.04"
            )

        elif prebuilt_genome == GenomeType.drosophilia:
            reference_zipped_p = LatchFile(
                "s3://latch-public/test-data/18440/pipseeker-gex-reference-dm-flybase-r6-v47-2022.09.tar.gz"
            ).local_path
            reference_p = Path(
                "/root/pipseeker-gex-reference-dm-flybase-r6-v47-2022.09"
            )

        elif prebuilt_genome == GenomeType.zebrafish:
            reference_zipped_p = LatchFile(
                "s3://latch-public/test-data/18440/zebrafish_danio_rerio_GRCz11_r110_2023.08.tar.gz"
            ).local_path
            reference_p = Path("/root/zebrafish_danio_rerio_GRCz11_r110_2023.08")

        elif prebuilt_genome == GenomeType.arabidopsis_thaliana:
            reference_zipped_p = LatchFile(
                "s3://latch-public/test-data/18440/pipseeker-gex-reference-arabidopsis-thaliana-TAIR10.55-protein-coding-2023.02.tar.gz"
            ).local_path
            reference_p = Path(
                "/root/pipseeker-gex-reference-arabidopsis-thaliana-TAIR10.55-protein-coding-2023.02"
            )

        if get_path_only:
            return reference_zipped_p

        else:
            # Unpack the prebuilt reference.
            subprocess.run(
                ["tar", "-zxvf", f"{reference_zipped_p}", "-C", "/root"],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )

    # Using a custom prebuilt mapping reference.
    elif genome_source == "custom_prebuilt_genome":
        if custom_prebuilt_genome is not None:
            reference_p = Path(custom_prebuilt_genome)
            if get_path_only:
                return reference_p

        elif custom_prebuilt_genome_zipped is not None:
            reference_zipped_p = Path(custom_prebuilt_genome_zipped)
            if get_path_only:
                return reference_zipped_p

            reference_p = Path(f"/root/{reference_zipped_p.stem}").with_suffix("")

            print("Unpacking the custom prebuilt genome")
            unpacked_data = False  # Tracks whether the untar/unzip operation was attempted.
            if reference_zipped_p.suffixes[-2:] == [".tar", ".gz"]:
                subprocess.run(
                    ["tar", "-zxvf", str(reference_zipped_p), "-C", "/root"],
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
                unpacked_data = True

            elif reference_zipped_p.suffix == ".zip":
                subprocess.run(
                    ["unzip", "-o", str(reference_zipped_p), "-d", "/root"],
                    check=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
                unpacked_data = True

            # Check whether the uncompressed directory exists and matches the expected name.
            if unpacked_data:
                #   This handles the case where users might not pack their data inside of a single directory
                #   or the top-level zipped dir might use a different name than the parent zipped filename.
                #   If the directory is not found, the user will be prompted to re-upload the data.
                # List the files in the root dir and see if the directory name is present from the zipped file output
                root_dir = Path("/root")
                root_dir_files = [f.name for f in root_dir.iterdir()]

                # Check if the directory name is present in the root dir
                if not reference_zipped_p.stem in root_dir_files:
                    raise ValueError(f"Unpacking failed. The directory {reference_p} was not found.\n"
                                     "Please ensure that you compressed your reference with a single top-level "
                                     "directory containing the reference genome. \n"
                                     "Also ensure the top-level folder matches the prefix of your compressed file.")

            else:
                # Data was not unpacked, due to file extension mismatch.
                print('The provided genome must be compressed using .tar.gz or .zip format '
                      'or uploaded without compression.')
    else:
        print("No reference genome provided. Continuing.")
        print(f"Genome source: {genome_source}")
        reference_p = None

    return reference_p