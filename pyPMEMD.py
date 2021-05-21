#!/usr/bin/python
'''

Python Script For AMBER Molecular Dynamics Simulation
主程序 Version 1.03

Author YH. W
Last Update: 2021/05/20

'''

import os
import re
import sys
cwd = str(os.getcwd())
sys.path.append(cwd)

def get_pdb():
    '''
    尝试自动获取PDB ID 否则提示用户手动输入
    
    Return
    ----------
    PDB ID字符串
    '''

    def check(pdb):
        return re.fullmatch(r'^\d[0-9a-zA-Z]{3,}$', pdb)

    pdb = str(os.path.split(cwd)[1])
    if check(pdb):
        return pdb
    else:
        while True:
            pdb = str(
                input('\n要自动获取PDB代号 请将晶体文件夹修改为晶体ID\n请手动输入晶体PDB代号:')).strip().upper()
            if check(pdb):
                return pdb
            else:
                print('请输入正确的PDB代号!')


def get_wat_num(pdb):

    '''
    抓取主链与非主链分界Residue Number
    
    Paramters
    ----------
    pdb: PDB ID字符串

    Return
    ----------
    分界处（主链 | 非主链） 非主链Residue Number
    '''

    return os.popen("cat %scomsolvate.pdb | grep WAT -w | awk '{print $5}'" % pdb).readlines()


def mina_prepare(wat_start):
    '''
    约束主链能量最小化 输入文件编写
    
    Paramters
    ----------
    wat_start: 非主链起始Residue Number    
    '''

    text = ''
    text += '&cntrl\n'
    text += 'imin=1,\n' #能量最小化模式
    text += 'cut=8.0,\n' #非键截短值 8埃
    text += 'ntpr=50,\n'    #每50步记录一次信息
    text += 'ntb=1,\n'  
    text += 'ntr=1,\n'  #允许约束原子
    text += 'maxcyc=10000,\n'   #最小化的最大循环数
    text += 'ncyc=5000,\n'  #在5000步后 优化方法由最陡下降法转变为共轭梯度法
    text += 'restraint_wt=2.0,\n'
    text += 'restraintmask='    #约束的原子标识
    text += "':1-%s'\n" % str(wat_start-1)  #约束主链
    text += '/'
    text += '\nEND'
    with open('./mina/mina.in', 'w') as f:
        f.write(text)


def minb_prepare(wat_start, wat_end):
    '''
    约束非主链能量最小化 输入文件编写
    
    Paramters
    ----------
    wat_start: 非主链起始Residue Number
    wat_end: 非主链末尾Residue Number
    '''

    text = ''
    text += '&cntrl\n'
    text += 'imin=1,\n'
    text += 'cut=8.0,\n'
    text += 'ntpr=50,\n'
    text += 'ntb=1,\n'
    text += 'ntr=1,\n'
    text += 'maxcyc=10000,\n'
    text += 'ncyc=5000,\n'
    text += 'restraint_wt=2.0,\n'
    text += 'restraintmask='
    text += "':%s-%s'\n" % (str(wat_start), str(wat_end))   #约束水分子
    text += '/'
    text += '\nEND'
    with open('./minb/minb.in', 'w') as f:
        f.write(text)


def minc_prepare():
    '''
    无约束能量最小化 输入文件编写
    '''

    text = ''
    text += '&cntrl\n'
    text += 'imin=1,\n'
    text += 'cut=8.0,\n'
    text += 'ntpr=50,\n'
    text += 'ntb=1,\n'
    text += 'maxcyc=10000,\n'
    text += 'ncyc=5000,\n'
    text += '/' #无约束最小化
    text += '\nEND'
    with open('./minc/minc.in', 'w') as f:
        f.write(text)


def nvt_prepare():
    '''
    体系恒容恒压加热100ps 输入文件编写
    '''

    text = ''
    text += '&cntrl\n'
    text += 'imin=0,\n' #分子动力学模拟模式
    text += 'tempi=0.0, temp0=300, ntt=3, gamma_ln=2.0,\n'  #初始温度0K 维持温度300K 使用Langevin算法 伽马碰撞频率2.0/ps
    text += 'ntp=1, taup=2.0,\n'    #恒压模式 压力松弛时间 2ps
    text += 'ntb=2,ntc=2,ntf=2,\n'      #ntb恒压模式 ntc约束与氢原子相关的键 ntf忽略氢键作用
    text += 'nstlim=50000, dt=0.002,\n' #MD步数50000步 步长0.002ps 总模拟时间=50000*0.002=100ps
    text += 'ntx=1,irest=0,cut=8.0,\n'  #从inpcrd读取坐标等参数 启动新模拟而非重启模拟
    text += 'ntpr=1000, ntwr=1000, ntwx=1000,\n'    #每1000步记录能量 坐标 重启文件信息
    text += '/'
    text += "&wt TYPE='TEMP0',ISTEP1=0, ISTEP2=25000,\n"
    text += 'VALUE1=0.0, VALUE2=300.0, /\n'
    text += "&wt TYPE='END' /\n"
    text += 'END'
    with open('./nvt/nvt.in', 'w') as f:
        f.write(text)


def npt_prepare():
    '''
    Molecular Dynamics主模拟 输入文件编写
    '''
    
    text = ''
    text += '&cntrl\n'
    text += 'imin=0,\n'
    text += 'temp0=300, ntt=3, gamma_ln=5.0,\n' #伽马碰撞频率有变动
    text += 'ntp=1, taup=2.0,\n'
    text += 'ntb=2,ntc=2,ntf=2,\n'
    text += 'nstlim=50000000, dt=0.002,\n'  #MD步数50000000步 步长0.002ps 总模拟时间=50000000*0.002=100000ps=100ns
    text += 'ntx=5,irest=1,\n'  #继续加热体系计算
    text += 'ntpr=1000, ntwr=1000, ntwx=1000,\n'
    text += '/'
    text += '\nEND'
    with open('./npt/npt.in', 'w') as f:
        f.write(text)

def mkdir(path):
    '''
    创建必要文件夹

    Paramters
    ----------
    path: 文件夹PATH
    '''

    try:
        os.makedirs(path)
    except FileExistsError:
        pass

def error_check(code):
    '''
    检查是否发生错误并中断MD流程

    Paramters
    ----------
    code: shell命令执行返回码
    '''

    if code != 0:
        print('Error: Error Occurred When Current Task Runing!\n')
        sys.exit(2)

def main():
    '''
    启动Molecular Dynamics Simulation 主函数
    '''

    pdb = get_pdb()

    mkdir('./mina')
    mkdir('./minb')
    mkdir('./minc')
    mkdir('./nvt')
    mkdir('./npt')

    wat_num = get_wat_num(pdb)
    wat_start = int(wat_num[0])
    wat_end = int(wat_num[-1])

    mina_prepare(wat_start)
    minb_prepare(wat_start, wat_end)
    minc_prepare()
    nvt_prepare()
    npt_prepare()
    print('\nStart to MD Process...')

    mina_cmd = 'pmemd.cuda -O -i ./mina/mina.in -o ./mina/mina.out -p %scomsolvate.prmtop -c %scomsolvate.inpcrd -r ./mina/mina.rst -ref %scomsolvate.inpcrd' % (pdb,pdb,pdb)
    minb_cmd = 'pmemd.cuda -O -i ./minb/minb.in -o ./minb/minb.out -p %scomsolvate.prmtop -c ./mina/mina.rst -r ./minb/minb.rst -ref ./mina/mina.rst' % pdb
    minc_cmd = 'pmemd.cuda -O -i ./minc/minc.in -o ./minc/minc.out -p %scomsolvate.prmtop -c ./minb/minb.rst -r ./minc/minc.rst' % pdb
    nvt_cmd = 'pmemd.cuda -O -i ./nvt/nvt.in -o ./nvt/nvt.out -p %scomsolvate.prmtop -c ./minc/minc.rst -r ./nvt/nvt.rst -x ./nvt/nvt.mdcrd' % pdb
    npt_cmd = 'pmemd.cuda -O -i ./npt/npt.in -o ./npt/npt.out -p %scomsolvate.prmtop -c ./nvt/nvt.rst -r ./npt/npt.rst -x ./npt/npt.mdcrd' % pdb

    print('\nRunning Minimize Progress A...\n')
    re1 = os.system(mina_cmd)
    error_check(re1)

    print('\nRunning Minimize Progress B...\n')
    re2 = os.system(minb_cmd)
    error_check(re2)

    print('\nRunning Minimize Progress C...\n')
    re3 = os.system(minc_cmd)
    error_check(re3)

    print('\nRunning 100ps Heating Progress...\n')
    re4 = os.system(nvt_cmd)
    error_check(re4)

    print('\nRunning 100ns Molecular Dynamics Simulation...\n')
    re5 = os.system(npt_cmd)
    error_check(re5)

    print('\nMolecular Dynamic Simulation Complete.')
    print(''.center(80,'-'))


if __name__ == '__main__':
    main()
