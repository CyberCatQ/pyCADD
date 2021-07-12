#!/usr/bin/python
'''

MD轨迹最低能量（势能）构象提取 
Version 1.10

Author YH. W
Last update: 2021/07/12

'''

import os
import re
import sys
import pytraj as pt
root_path = os.path.abspath(os.path.dirname(__file__)).split('src')[0]  # 项目路径 绝对路径
pdb_path = root_path.split('automatedMD')[0]                            # PDB项目绝对路径(如果有)
libpynalysis_path = root_path + '/src/libpynalysis/'                    # libpynalysis文件夹路径

try:
    os.makedirs(pdb_path + 'pynalysis/lowestenergy/')
except FileExistsError:
    pass

if not os.path.exists(libpynalysis_path + 'lowestenergy/process_mdout.perl'):
    raise RuntimeError('\nprocess_mdout.perl is missing.\nPlease place it in automatedMD/src/libpynalisis/lowestenergy and try again.')
else:
    os.system('cp %s/lowestenergy/process_mdout.perl %spynalysis/lowestenergy' % (libpynalysis_path, pdb_path))

def main(traj):
    '''
    最低能量（势能）构象帧提取

    Parameters
    ----------
    traj : <object pytraj.traj>
        Pytraj MD轨迹对象
    '''

    print('\nStart Lowest Energy Structure Analysis...')

    # 调用perl脚本分析MD out数据
    # process_mdout.perl脚本为AMBER官方脚本
    os.system(
        'cd %spynalysis/lowestenergy/ && chmod 777 ./process_mdout.perl && ./process_mdout.perl %snpt/npt.out' % (pdb_path, pdb_path))
    print('\nProcessing Complete.', end='\n')

    print('\nSearching Lowest Energy Time...', end='\n')
    temp = os.popen('''
    cat  %spynalysis/lowestenergy/summary.EPTOT | awk '{if($2<min) {min=$2;print $1"   "min}}'
    ''' % pdb_path)    #EPTOT总势能

    lines = temp.readlines()
    match_time = re.match('\d{1,5}.000', lines[-1].strip())  # 匹配最后一行（最低能量）的时间
    lowest_energy_time = match_time.group()  # 获取时间具体值

    for info in lines:
        print(info.strip())

    print('\nLowest Energy Time:', lowest_energy_time, 'ps')

    os.environ['lowest_energy_time'] = str(
        lowest_energy_time)  # 向shell传递变量 搜索最低能量帧
    temp2 = os.popen('grep $lowest_energy_time %snpt/npt.out' % pdb_path)
    frame_temp = temp2.read().strip()  # 最低能量帧信息
    print('\nLowest Energy Frame info: \n', frame_temp)

    # 提取最低能量帧数: 帧数 = 步数 / 1000
    match_frame = re.search('(?<=NSTEP\s=\s)\d*(?=000)', str(frame_temp))
    lowest_energy_frame = int(match_frame.group())
    print('\nLowest Energy Frame:', lowest_energy_frame)

    pt.write_traj('%spynalysis/lowestenergy/lowest_energy_structure.pdb' % pdb_path, traj,
                  frame_indices=[lowest_energy_frame], overwrite=True)  # 输出最低能量构象PDB文件
    print('\nLowest Energy Structure PDB file saved.')


if __name__ == '__main__':

    print('请从主程序pynalysis.py中运行此脚本')
