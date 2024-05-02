import unittest
import os
from latch.types import LatchDir, LatchFile
from wf.resource_estimator import get_num_threads, get_memory_requirement_gb, get_disk_requirement_gb, mapping_ref_size_estimator
from wf.configurations import GenomeType
from unit_tests.test_utils import UnitTest


class ResourceEstimatorTest(UnitTest):
    def setUp(self):
        super().setUp()
        # Check whether the fastqs exist in latch account.
        self.latch_fastq_dir_path = 'latch://26230.account/TEST_fastqs_10mb_Pbmc_v4'
        self.fastqs_local_path = os.path.join(self.TEST_DATA, 'fastqs')

        # Prep test fastqs on the Latch cloud.
        try:
            LatchDir(self.latch_fastq_dir_path)
            print(f"Successfully detected test FASTQs on Latch at {self.latch_fastq_dir_path}")
        except:
            # Upload files to latch via bash.
            print(f'Uploading test fastqs to Latch at {self.latch_fastq_dir_path}')
            os.system(f'latch cp {self.fastqs_local_path} {self.latch_fastq_dir_path}')
            if not os.path.exists(LatchDir(self.latch_fastq_dir_path)):
                raise FileNotFoundError('Failed to upload or identify fastqs to Latch.')

        # Prep STAR index on the Latch cloud.
        # Zipped STAR index
        self.star_index_zipped_path_local = os.path.join(self.TEST_DATA, 'STAR_test_index.zip')
        self.latch_star_zipped_path = 'latch://26230.account/STAR_test_index.zip'

        try:
            LatchFile(self.latch_star_zipped_path)
            print(f"Successfully detected test STAR index on Latch at {self.latch_star_zipped_path}")
        except:
            print(f'Uploading STAR index to Latch at {self.latch_star_zipped_path}')
            os.system(f'latch cp {self.star_index_zipped_path_local} {self.latch_star_zipped_path}')
            if not os.path.exists(LatchFile(self.latch_star_zipped_path)):
                raise FileNotFoundError('Failed to upload or identify STAR index to Latch.')

        # Unzipped STAR index
        self.star_index_dir_path_local = os.path.join(self.TEST_DATA, 'STAR_test_index')
        self.latch_star_unzipped_dir_path = 'latch://26230.account/STAR_test_index'

        try:
            LatchDir(self.latch_star_unzipped_dir_path)
            print(f"Successfully detected test STAR index on Latch at {self.latch_star_unzipped_dir_path}")
        except:
            #  First, unzip the STAR index locally.
            print(f'Unzipping STAR index to {self.star_index_dir_path_local}')
            os.system(f'unzip {self.star_index_zipped_path_local} -d {self.star_index_dir_path_local}')
            if not os.path.exists(self.star_index_dir_path_local):
                raise FileNotFoundError('Failed to unzip STAR index.')
            print(f'Uploading STAR index to Latch at {self.latch_star_unzipped_dir_path}')
            os.system(f'latch cp {self.star_index_dir_path_local} {self.latch_star_unzipped_dir_path}')
            if not os.path.exists(LatchDir(self.latch_star_unzipped_dir_path)):
                raise FileNotFoundError('Failed to upload or identify STAR index to Latch.')

    def test_get_num_threads(self):
        # Dir has very tiny fastq set so should have the 4-thread minimum.
        self.assertEqual(get_num_threads(fastq_directory=LatchDir(self.latch_fastq_dir_path)), 4)

    def test_get_mapping_ref_size(self):

        # Zipped STAR index.
        ref_size_bytes = mapping_ref_size_estimator(genome_source='custom_prebuilt_genome',
                                                prebuilt_genome=None,
                                                custom_prebuilt_genome=None,
                                                custom_prebuilt_genome_zipped=LatchFile(self.latch_star_zipped_path))
        self.assertEqual(ref_size_bytes, 1965038.75)  # 0.00183 GB

        # Unzipped STAR index.
        ref_size_bytes = mapping_ref_size_estimator(genome_source='custom_prebuilt_genome',
                                                prebuilt_genome=None,
                                                custom_prebuilt_genome=LatchDir(self.latch_star_unzipped_dir_path),
                                                custom_prebuilt_genome_zipped=None)
        self.assertEqual(ref_size_bytes, 2169065)  # 0.00202 GB

        # For the pre-built references hosted on s3.
        #   Human ref is 9.61 GB.
        ref_size = mapping_ref_size_estimator(genome_source='prebuilt_genome',
                                              prebuilt_genome=GenomeType.human,
                                              custom_prebuilt_genome=None,
                                              custom_prebuilt_genome_zipped=None)
        ref_size_gb = ref_size / 1024 ** 3
        self.assertAlmostEqual(ref_size_gb, 9.61, places=2)

    def test_get_memory_requirement_gb(self):
        # Zipped STAR index.
        required_ram = get_memory_requirement_gb(fastq_directory=LatchDir(self.latch_fastq_dir_path),
                                                 genome_source='custom_prebuilt_genome',
                                                 prebuilt_genome=None,
                                                 custom_prebuilt_genome=None,
                                                 custom_prebuilt_genome_zipped=LatchFile(self.latch_star_zipped_path),
                                                 downsample_to=None,
                                                 input_reads=None,
                                                 sorted_bam=False)
        self.assertEqual(required_ram, 22)

        # Unzipped STAR index.
        required_ram = get_memory_requirement_gb(fastq_directory=LatchDir(self.latch_fastq_dir_path),
                                                 genome_source='custom_prebuilt_genome',
                                                 prebuilt_genome=None,
                                                 custom_prebuilt_genome=None,
                                                 custom_prebuilt_genome_zipped=LatchDir(
                                                     self.latch_star_unzipped_dir_path),
                                                 downsample_to=None,
                                                 input_reads=None,
                                                 sorted_bam=False)
        self.assertEqual(required_ram, 22)

        # For the pre-built references hosted on s3.
        #   Human ref is 9.61 GB.
        required_ram = get_memory_requirement_gb(fastq_directory=LatchDir(self.latch_fastq_dir_path),
                                               snt_fastq=None,
                                               hto_fastq=None,
                                               genome_source='prebuilt_genome',
                                               prebuilt_genome=GenomeType.human,
                                               custom_prebuilt_genome=None,
                                               custom_prebuilt_genome_zipped=None,
                                               downsample_to=None,
                                               input_reads=None,
                                               sorted_bam=False)
        self.assertEqual(required_ram, 22)

    def test_disk_requirement_gb(self):
        # On s3 (human genome is 9.61 GB).
        required_ram = get_disk_requirement_gb(fastq_directory=LatchDir(self.latch_fastq_dir_path),
                                               snt_fastq=None,
                                               hto_fastq=None,
                                               genome_source='prebuilt_genome',
                                               prebuilt_genome=GenomeType.human,
                                               custom_prebuilt_genome=None,
                                               custom_prebuilt_genome_zipped=None,
                                               downsample_to=None,
                                               input_reads=None,
                                               sorted_bam=False)
        self.assertEqual(required_ram, 32)

        # Zipped STAR index.
        #   Will use 2GB for the zipped index, since STAR index is only 0.002 GB and rounds to 0.
        required_disk = get_disk_requirement_gb(fastq_directory=LatchDir(self.latch_fastq_dir_path),
                                                genome_source='custom_prebuilt_genome',
                                                prebuilt_genome=None,
                                                custom_prebuilt_genome=None,
                                                custom_prebuilt_genome_zipped=LatchFile(self.latch_star_zipped_path),
                                                downsample_to=None,
                                                input_reads=None,
                                                sorted_bam=False)
        self.assertEqual(required_disk, 2)

        # Unzipped STAR index.
        #   Again, uses 2GB because of small unzipped index size.
        required_disk = get_disk_requirement_gb(fastq_directory=LatchDir(self.latch_fastq_dir_path),
                                                genome_source='custom_prebuilt_genome',
                                                prebuilt_genome=None,
                                                custom_prebuilt_genome=None,
                                                custom_prebuilt_genome_zipped=LatchDir(
                                                    self.latch_star_unzipped_dir_path),
                                                downsample_to=None,
                                                input_reads=None,
                                                sorted_bam=False)
        self.assertEqual(required_disk, 2)

        # No fastqs.
        required_disk = get_disk_requirement_gb(fastq_directory=None,
                                                genome_source='custom_prebuilt_genome',
                                                prebuilt_genome=None,
                                                custom_prebuilt_genome=None,
                                                custom_prebuilt_genome_zipped=LatchDir(
                                                    self.latch_star_unzipped_dir_path),
                                                downsample_to=None,
                                                input_reads=None,
                                                sorted_bam=False)
        self.assertEqual(required_disk, 2)

        # No fastqs, no genome.
        required_disk = get_disk_requirement_gb(fastq_directory=None,
                                                genome_source=None,
                                                prebuilt_genome=None,
                                                custom_prebuilt_genome=None,
                                                custom_prebuilt_genome_zipped=None,
                                                downsample_to=None,
                                                input_reads=None,
                                                sorted_bam=False)
        self.assertEqual(required_disk, 2)


if __name__ == '__main__':
    unittest.main()
