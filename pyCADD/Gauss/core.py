import logging
import os
from rich.prompt import Confirm

from pyCADD.utils.check import check_file
from pyCADD.utils.tool import tail_progress

logger = logging.getLogger(__name__)

def generate_opt(original_st:str, charge:int, multiplicity:int, dft:str='B3LYP', basis_set:str='6-31g*', solvent:str='water', loose:bool=True, correct:bool=True, td:bool=False):
    '''
    生成Gaussian结构优化输入文件
    
    Parameters
    ----------
    original_st : str
        原始分子结构文件路径
    charge : int 
        电荷量
    multiplicity : int 
        自旋多重度
    dft : str
        泛函数
    basis_set : str
        基组
    solvent : str
        PCM模型溶剂
    loose : bool
        是否提高优化任务中的收敛限 更快收敛
    td : bool
        是否为激发态结构优化计算(计算荧光/磷光发射能用)

    Return
    ----------
    str
        生成的输入文件名称
    '''
    molname = os.path.basename(os.path.abspath(original_st)).split('.')[0]
    if td:
        td_suffix = '_td'
    else:
        td_suffix = ''

    opt_file = molname + '_opt%s.gjf' % td_suffix
    chk_file = molname + '_opt%s.chk' % td_suffix

    if check_file(opt_file):
        if not Confirm.ask('%s is existed. Overwrite?' % opt_file, default=True):
            return opt_file, chk_file

    # 调用Multiwfn预生成输入文件
    os.system('''
    Multiwfn %s > /dev/null << EOF
    100
    2
    10
    tmp.gjf
    0
    q
    EOF
    ''' % original_st)
    # 溶剂可为None
    if solvent == 'None':
        scrf = ''
    else:
        scrf = ' scrf(solvent=%s)' % solvent
    # 提高收敛限
    if loose:
        opt_config = 'opt=loose'
    else:
        opt_config = 'opt'
    # 激发态
    if td:
        TD = ' TD'
    else:
        TD = ''
    # 色散矫正
    if correct:
        correct_cofig = ' em=GD3BJ'
    else:
        correct_cofig = ''

    keyword = '# %s %s/%s%s%s%s' % (opt_config, dft, basis_set, correct_cofig, TD, scrf)

    os.system('''
    cat << EOF > %s
%s=%s
%s

%s optimize

  %s %s
EOF
''' % (opt_file, r"%chk", chk_file, keyword, molname + td_suffix, charge, multiplicity))

    os.system("awk '{if (NR>5) print }' tmp.gjf >> %s" % opt_file)
    os.remove('tmp.gjf')

    return opt_file, chk_file

def generate_energy(original_st:str, charge:int, multiplicity:int, dft:str='B3LYP', basis_set:str='6-31g*', solvent:str='water', correct:bool=True, td:bool=False):
    '''
    生成Gaussian单点能计算输入文件
    
    Parameters
    ----------
    original_st : str
        原始分子结构文件路径
    charge : int 
        电荷量
    multiplicity : int 
        自旋多重度
    dft : str
        泛函数
    basis_set : str
        基组
    solvent : str
        PCM模型溶剂

    Return
    ----------
    str
        生成的输入文件名称
    '''
    molname = os.path.basename(os.path.abspath(original_st)).split('.')[0]
    if td:
        td_suffix = '_td'
    else:
        td_suffix = ''

    energy_file = molname + '_energy%s.gjf' % td_suffix
    chk_file = molname + '_energy%s.chk' % td_suffix
    if check_file(energy_file):
        if not Confirm.ask('%s is existed. Overwrite?' % energy_file, default=True):
            return energy_file, chk_file

    # 调用Multiwfn预生成输入文件
    os.system('''
    Multiwfn %s > /dev/null << EOF
    100
    2
    10
    tmp.gjf
    0
    q
    EOF
    ''' % original_st)

    # 溶剂可为None
    if solvent == 'None':
        scrf = ''
    else:
        scrf = ' scrf(solvent=%s)' % solvent
    # 激发态
    if td:
        TD = ' TD'
    else:
        TD = ''
    # 色散矫正
    if correct:
        correct_cofig = ' em=GD3BJ'
    else:
        correct_cofig = ''

    keyword = '# %s/%s%s%s%s' % (dft, basis_set, correct_cofig, TD, scrf)

    os.system('''
    cat << EOF > %s
%s=%s
%s

%s Single Point Energy

  %s %s
EOF
''' % (energy_file, r"%chk", chk_file, keyword, molname + td_suffix, charge, multiplicity))

    os.system("awk '{if (NR>5) print }' tmp.gjf >> %s" % energy_file)
    os.remove('tmp.gjf')

    return energy_file, chk_file

def get_gaussian():
    '''
    获取高斯可执行文件路径
    Return
    ---------
    str
        高斯可执行文件路径
    '''

    g16 = os.popen('which g16').read().strip()
    g09 = os.popen('which g09').read().strip()
    if g16:
        gaussian = g16
    elif g09:
        logger.warning('You are using gaussian 09, that may cause some unknown errors.\nGaussian 16 is recommend.')
        gaussian = g09
    else:
        logger.error('Gaussian is not installed.')
        return None
    return gaussian
    

def system_default(gauss_path:str, cpu_count:int, memory:str):
    '''
    修改系统计算资源占用设定

    Parameters
    ----------
    gauss_path : str
        高斯可执行文件路径
    cpu_count : int 
        CPU核心使用数量
    memory : str
        内存占用大小(MB/GB)

    Return
    ----------
    str
        Gaussian 可执行文件路径
    '''

    default_route = os.path.dirname(gauss_path) + '/Default.Route'
    with open(default_route,'w') as f:
        f.write('-P- %s\n' % cpu_count)
        f.write('-M- %s\n' % memory)
    logger.debug('Default system setting changed: CPU = %s  Mem = %s' % (cpu_count, memory))
    
def generate_fchk(chk_file:str):
    '''
    生成fchk文件
    '''
    os.system('formchk %s' % chk_file)
    logger.debug('fchk file %s is saved.' % chk_file.split('.')[0] + '.fchk')


def _check_gauss_finished(line):
    '''
    检查高斯计算任务是否结束

    Parameter
    ---------
    line : str
        高斯输出文件单行内容

    Return
    ---------
    bool
        任务结束返回True 否则返回False
    '''

    if 'Normal termination of Gaussian' in line:
        return True
    else:
        return False

def tail_gauss_job(output_file:str):
    '''
    追踪高斯计算任务进度
    Parameter
    ---------
    output_file : str
        高斯计算任务输出文件路径
    '''
    tail_progress(output_file, _check_gauss_finished)
    

