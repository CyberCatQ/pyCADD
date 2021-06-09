import csv
import getopt
import multiprocessing
import os
import re
import sys
from getopt import getopt

from schrodinger import structure as struc
from schrodinger.application.glide import poseviewconvert as pvc
from schrodinger.job import jobcontrol as jc
from schrodinger.protein import getpdb

cwd = str(os.getcwd())
sys.path.append(cwd)


def error_handler(error):
    print(error.__cause__)
    print(''.center(80, '-'), end='\n')


def success_handler(result):
    global now_complete
    now_complete += result
    print("%s Job(s) Successd / Total: %s" % (now_complete, total))


class Console:

    def __init__(self, pdbid=None, pdbfile=None, ligname=None, flag=None) -> None:
        self.pdbid = pdbid
        self.pdbfile = pdbfile
        self.ligname = ligname
        self.flag = flag

    @staticmethod
    def checkpdb(self,pdb):
        '''
        检查PDB ID合法性

        Return
        ----------
        合法返回匹配对象 否则返回None
        '''
        return re.fullmatch(r'^\d[0-9a-zA-Z]{3,}$', pdb)

    def get_pdbid(self):
        '''
        获取用户输入的PDBID并检查合法性

        Return
        ----------
        PDB ID字符串

        '''

        if not self.pdbid:
            pdb = str(os.path.split(cwd)[1])

        if self.checkpdb(pdb):
            self.pdbid = pdb
            return pdb
        else:
            while True:
                pdb = str(
                    input('\n要自动获取PDB ID 请将晶体文件夹修改为PDB ID\n请手动输入晶体PDB ID:')).strip().upper()
                if self.checkpdb(pdb):
                    self.pdbid = pdb
                    return pdb
                else:
                    print('请输入正确的PDB ID!')

    def _preprocess(self):
        '''
        一些预处理操作：
        下载PDB文件 | 检查PDB晶体类型 Apo/单体/多聚体 | 是否保留单链

        Returen
        ----------
        选择保留单链得到的文件PATH
        '''
        if not self.pdbfile:
            self.pdbfile = self.pdbid + '.pdb'

        if not os.path.exists(self.pdbfile):  # 下载PDB文件
            getpdb.get_pdb(self.pdb)

        lig_lis = os.popen(
            "cat %s.pdb | grep -w -E ^HET | awk '{print $2}'" % pdb).readlines()

        if len(lig_lis) == 0:
            print('%s为Apo蛋白晶体 无法自动处理.' % pdb)
            sys.exit(1)
        elif len(lig_lis) > 1:
            print('\n')
            os.system('cat %s.pdb | grep -w -E ^HET' % pdb)
            print('存在多个配体小分子 是否需要保留单链？(Y/N)')
            k = input().strip().upper()

            if k == 'Y':
                file = pdb + '.pdb'
                chain = input('输入保留链名:').strip().upper()
                pdb_file = keep_chain(pdb, file, chain)
                return pdb_file
            elif k == 'N':
                print('\nChoosed Intact Crystal.\n')
                return pdb_file
        else:
            return pdb_file


def main():
    console = Console()
    console.get_pdbid()