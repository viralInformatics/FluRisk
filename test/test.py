# !/usr/bin/python
# -*- coding: utf-8 -*-
# Author:lihuiru
# Created on 2024/3/13 11:01
import os
import subprocess
import sys
import traceback
import unittest
from contextlib import redirect_stdout, redirect_stderr

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class FlupreTest(unittest.TestCase):
    directory_path = "./"

    @classmethod
    def setUpClass(cls):
        # 在所有测试开始之前执行一次，用于设置测试所需的全局配置

        cls.test_files_path = os.path.join(cls.directory_path, 'test_files')
        cls.standardized_fasta_path = os.path.join(cls.directory_path, 'standardized_fasta')
        cls.adaptation_path = os.path.join(cls.directory_path, 'adaptation')
        cls.virulence_path = os.path.join(cls.directory_path, 'virulence')
        cls.binding_path = os.path.join(cls.directory_path, 'binding')
        cls.result_path = os.path.join(cls.directory_path, 'result')
        cls.ada_prediction_path = os.path.join(cls.directory_path, 'ada_prediction')
        cls.vir_prediction_path = os.path.join(cls.directory_path, 'vir_prediction')
        cls.bin_prediction_path = os.path.join(cls.directory_path, 'bin_prediction')

        # 创建所有需要的目录
        os.makedirs(cls.result_path, exist_ok = True)
        os.makedirs(cls.ada_prediction_path, exist_ok = True)
        os.makedirs(cls.vir_prediction_path, exist_ok = True)

    def run_command(self, command):
        # 辅助函数，用于执行命令并返回结果 (modified to support Python 3.6)
        process = subprocess.run(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE,
                                 universal_newlines = True)
        return process

    def compress_marker_directories(self):
        """
        压缩所有marker类型的输出目录为zip文件
        """
        marker_types = ['adaptation', 'binding', 'resistance', 'virulence', 'transmissibility']

        for marker_type in marker_types:
            if os.path.exists(marker_type) and os.listdir(marker_type):
                try:
                    # 使用系统zip命令，递归压缩目录
                    compress_result = subprocess.run(['zip', '-r', f'{marker_type}/{marker_type}.zip', marker_type],
                                                     check = True, capture_output = True, text = True)
                    print(f"Created zip file: {marker_type}.zip")
                except subprocess.CalledProcessError as e:
                    print(f"Error creating zip for {marker_type}: {e}")
                except FileNotFoundError:
                    print(f"zip command not found, skipping compression for {marker_type}")
            else:
                print(f"Directory {marker_type} does not exist or is empty, skipping zip creation.")

    def test_anno_command(self):
        # 测试anno功能
        command = f'flupre anno -i {self.test_files_path} -o {self.result_path}'
        self.test_output = self.run_command(command)  # 修改这里，保存输出到self.test_output
        self.assertEqual(self.test_output.returncode, 0, "anno command failed with error: " + self.test_output.stderr)

    def test_extract_command(self):
        # 测试extract功能
        command = f'flupre extract -i {self.standardized_fasta_path} -a {self.result_path}'
        self.test_output = self.run_command(command)
        self.assertEqual(self.test_output.returncode, 0,
                         "extract command failed with error: " + self.test_output.stderr)

        # # extract完成后立即压缩marker目录
        # print("Compressing marker directories after extract...")
        # self.compress_marker_directories()

    def test_predh_command(self):
        # 测试predh功能
        command = f'flupre predh -i {self.adaptation_path} -o {self.ada_prediction_path}'
        self.test_output = self.run_command(command)
        self.assertEqual(self.test_output.returncode, 0, "predh command failed with error: " + self.test_output.stderr)

    def test_predv_command(self):
        # 测试predv功能
        command = f'flupre predv -i {self.virulence_path} -o {self.vir_prediction_path}'
        self.test_output = self.run_command(command)
        self.assertEqual(self.test_output.returncode, 0, "predv command failed with error: " + self.test_output.stderr)

    def test_predr_command(self):
        # 测试predr功能
        command = f'flupre predr -i {self.binding_path} -o {self.bin_prediction_path}'
        self.test_output = self.run_command(command)
        self.assertEqual(self.test_output.returncode, 0, "predr command failed with error: " + self.test_output.stderr)

    def test_to_bar_command(self): #加一个for是因为测试时调用函数是按照test后的字母顺序，要保证在extract和predh predv后调用测试
        command = f'barplot -i ./ -o total_barplot/'
        self.test_output = self.run_command(command)
        self.assertEqual(self.test_output.returncode, 0, "bar command failed with error: " + self.test_output.stderr)

    def test_to_radar_command(self):
        # command = f'radarplot -i ./ -o radar_plot/'
        command = f'risk2tab -i ./ -a {self.result_path}  -o radar_file/'
        self.test_output = self.run_command(command)
        self.assertEqual(self.test_output.returncode, 0, "radar command failed with error: " + self.test_output.stderr)

        # # extract完成后立即压缩marker目录
        # print("Compressing marker directories after extract...")
        # self.compress_marker_directories()

    def tearDown(self):
        # Method to perform cleanup after tests, including capturing logs and exceptions
        test_method_name = self.id().split('.')[-1]  # Getting the method name
        print(test_method_name)
        test_method_dir = os.path.join(self.directory_path, 'runlog', test_method_name)
        os.makedirs(test_method_dir, exist_ok = True)

        # If an exception occurred, capture the information
        exception_info = sys.exc_info()
        has_exception = exception_info[0] is not None

        log_path = os.path.join(test_method_dir, 'log.txt')

        # Capture the logs and exceptions (if any) in the log file
        with open(log_path, 'w') as log_file, redirect_stdout(log_file), redirect_stderr(log_file):
            if has_exception:
                traceback.print_exception(*exception_info, file = log_file)

            # Capture the stdout and stderr from the test method (if any)
            test_output = getattr(self, 'test_output', None)
            if test_output:
                stdout, stderr = test_output.stdout, test_output.stderr
                if stdout:
                    log_file.write("\nStandard Output:\n" + stdout)
                if stderr:
                    log_file.write("\nStandard Error:\n" + stderr)


if __name__ == '__main__':
    unittest.main()
