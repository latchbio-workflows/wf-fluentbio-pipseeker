"""
Usage:
    This script should be run individually after opening up a development environment as follows:
        # In command line:

        # Register the workflow
        cd /path/to/latch-project-containing-wf-directory
        latch register .

        # After registering is complete:
        latch develop .

        # Run the test.
        python3 unit_tests/test_pipseeker.py

    In the future, this test will be updated to run within the IDE.
"""

import unittest
import os
from wf.configurations import GenomeType, Chemistry
from wf import pipseeker_wf
from unit_tests.test_utils import UnitTest
import glob

class PIPseekerTest(UnitTest):

    def test_snt_mode(self):
        self.star_index_zipped_path_local = os.path.join(self.TEST_DATA, 'STAR_test_index.zip')
        self.local_snt_dir_path = os.path.join(self.TEST_DATA, 'snt_tests')

        # Test snt mode LOCAL
        pipseeker_wf(pipseeker_mode='full_mode',
                     chemistry=Chemistry.v4,
                     output_directory='.',
                     fastq_directory=os.path.join(self.local_snt_dir_path, 'fastqs'),
                     snt_fastq=os.path.join(self.local_snt_dir_path, 'snt_fastqs'),
                     snt_tags=os.path.join(self.local_snt_dir_path, 'tags_1_withExtra.csv'),
                     snt_position=3,
                     genome_source='custom_prebuilt_genome',
                     prebuilt_genome=GenomeType.human,
                     custom_prebuilt_genome=None,
                     custom_prebuilt_genome_zipped=self.star_index_zipped_path_local,
                     custom_genome_reference_gtf=None,
                     custom_genome_reference_fasta=None,
                     min_sensitivity=1,
                     max_sensitivity=3,
                     clustering_percent_genes=80,
                     diff_exp_genes=50)

        # Verify that the RNA analysis ran to completion, and that cell calling and clustering outputs
        # are based on RNA datre stored in the main output directory.
        sens_1_2_3 = ['sensitivity_1', 'sensitivity_2', 'sensitivity_3']
        self.output_dir = os.path.join('pipseeker_out')
        self.logs_dir = os.path.join(self.output_dir, 'logs')
        self.verify_log_file(included=['Saving molecule info file', 'Running clustering', 'sensitivity_3',
                                       'Creating summary report', 'Merging SNT data'],
                             excluded=['Merging HTO data'], msg='Failed')

        self.verify_dir_ls(os.path.join(self.output_dir, 'cell_calling'), sens_1_2_3, msg='Failed')
        self.verify_dir_ls(os.path.join(self.output_dir, 'clustering'), sens_1_2_3, msg='Failed')
        self.assertFalse(os.path.isdir(os.path.join(self.output_dir, 'SNT', 'cell_calling')), 'Failed')
        self.assertFalse(os.path.isdir(os.path.join(self.output_dir, 'SNT', 'clustering')), 'Failed')
        self.assertTrue(os.path.isfile(os.path.join(self.output_dir,'report.html')), 'Failed')

    def get_latest_log_file(self, logs_dir=None):
        if logs_dir is None:
            log_files = glob.glob(os.path.join(self.logs_dir, 'pipseeker*.log'))
        else:
            log_files = glob.glob(os.path.join(logs_dir, 'pipseeker*.log'))
        return sorted(log_files)[-1]

    def verify_log_file(self, logs_dir=None, **kwargs):
        log_file = self.get_latest_log_file(logs_dir=logs_dir)
        self.verify_file_content(log_file, **kwargs)

    def verify_file_content(self, filename, included=[], excluded=[], msg=None, exact=None, case_sensitive=True):
        # A helper to verify the content of any text file.
        # included and excluded can be a string or a list of strings.

        if type(included) is str:
            included = [included]
        if type(excluded) is str:
            excluded = [excluded]

        with open(filename) as fid:
            text = fid.read()

        if not case_sensitive:
            text = text.lower()

        for included_str in included:
            if not case_sensitive:
                included_str = included_str.lower()
            self.assertTrue(included_str in text,
                            msg='String "%s" not found in file %s, %s' % (included_str, filename, msg))
        for excluded_str in excluded:
            if not case_sensitive:
                excluded_str = excluded_str.lower()
            self.assertFalse(excluded_str in text,
                             msg='String "%s" unexpectedly found in file %s, %s' % (excluded_str, filename, msg))
        if exact is not None:
            self.assertEqual(text, exact, msg=msg)

    def verify_dir_ls(self, dir, expected, msg=None):
        # Verify directory contents against a list of expected files (names only, not full path).
        self.assertEqual(sorted(os.listdir(dir)), sorted(expected), msg=msg)

if __name__ == '__main__':
    unittest.main()
