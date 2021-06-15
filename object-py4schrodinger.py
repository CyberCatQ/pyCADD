import csv
import getopt
import multiprocessing
import os
import re
import sys
import getopt

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
    global total
    global now_complete
    now_complete += result
    print("%s Job(s) Successd / Total: %s" % (now_complete, total))


def load_st(st_file):
    '''
    读取结构

    Parameters
    ----------
    st_file: 需要读取结构的文件path 需要Schrodinger支持的文件格式

    Return
    ----------
    结构对象

    '''

    return next(struc.StructureReader(st_file))


class Console:

    def __init__(self, pdbid=None, ligname=None, flag=None) -> None:
        self.pdbid = pdbid
        self.ligname = ligname
        self.flag = flag
        self.pdbfile = self.pdbid + '.pdb'

    @staticmethod
    def launch(cmd):
        '''
        使用jobcontrol启动一项job并等待结束

        Parameters
        ----------
        cmd： 等待执行的命令字符串

        '''
        cmd_list = cmd.split(' ')  # launch_job以列表形式提交参数
        job = jc.launch_job(cmd_list)
        print('JobId: %s' % job.JobId, end='\n')
        print('Job Name: %s' % job.Name, end='\n')
        print('Job Status: %s' % job.Status, end='\n')
        job.wait()  # 阻塞进程 等待Job结束

    @staticmethod
    def checkpdb(pdb):
        '''
        检查PDB ID合法性

        Return
        ----------
        合法返回True 否则返回False
        '''
        match = re.fullmatch(r'^\d[0-9a-zA-Z]{3,}$', pdb)
        if match:
            return True
        else:
            return False

    @staticmethod
    def check_ligname(ligname): #可能没有必要
        '''
        检查配体名称合法性

        Return
        ----------
        合法返回True 否则返回False
        '''
        match = re.search('[0-9A-Z]{2,3}$', ligname)
        if match:
            return True
        else:
            return False

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

    def __keep_chain(self, chain_name):
        '''
        读取PDB晶体文件并将单一链的结构输出为pdb文件

        Parameters
        ----------
        pdb_code: PDB ID字符串
        origin_file: 待处理原始文件
        chain_name: 要保留的链名称

        Return
        ----------
        保留单链结构的文件PATH
        '''

        st = load_st(self.pdbfile)  # 读取原始PDB结构
        st_chain_only = st.chain[chain_name].extractStructure()
        file = '%s_chain_%s.mae' % (self.pdbid, chain_name)
        st_chain_only.write(file)
        return file

    def __preprocess(self):
        '''
        一些预处理操作：
        下载PDB文件 | 检查PDB晶体类型 Apo/单体/多聚体 | 是否保留单链

        Returen
        ----------
        选择保留单链得到的文件PATH
        '''

        pdbfile = self.pdbfile

        if not os.path.exists(pdbfile):  # 下载PDB文件
            getpdb.get_pdb(self.pdbid)

        lig_lis = os.popen(
            "cat %s | grep -w -E ^HET | awk '{print $2}'" % pdbfile).readlines()

        if len(lig_lis) == 0:
            raise RuntimeError('%s为Apo蛋白晶体 无法自动处理.' % self.pdbid)

        elif len(lig_lis) > 1:
            print('\n')
            os.system('cat %s.pdb | grep -w -E ^HET' % self.pdbid)
            print('存在多个配体小分子 是否需要保留单链？(Y/N)')
            _flag = input().strip().upper()

            if _flag == 'Y':
                chain = input('输入保留链名:').strip().upper()
                pdbfile = self.__keep_chain(chain)
                self.pdbfile = pdbfile  # [属性修改]修改pdbfile为单链文件
                return pdbfile
            else:
                print('\nChoosed Intact Crystal.\n')
                return pdbfile  # 主要晶体文件仍为原始晶体文件
        else:
            return pdbfile


def main():
    console = Console()
    console.get_pdbid()
