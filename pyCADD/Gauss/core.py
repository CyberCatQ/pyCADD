import atexit
import logging
import os
import re
import signal
import subprocess
import sys

from rich.prompt import Confirm

logger = logging.getLogger(__name__)
RATIO = 27.2114313131

def generate_opt(original_st: str, charge: int, multiplicity: int, dft: str = 'B3LYP', basis_set: str = '6-31g*', solvent: str = 'water', loose: bool = True, correct: bool = True, td: bool = False, freq:bool=False):
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
    freq : bool
        是否计算频率

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

    if os.path.exists(opt_file):
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
    # 频率
    if freq:
        freq_kw = ' freq'
    else:
        freq_kw = ''

    keyword = '# %s %s/%s%s%s%s%s' % (opt_config,
                                    dft, basis_set, correct_cofig, TD, scrf, freq_kw)

    os.system('''
    cat << EOF > %s
%s=%s
%s

%s optimize

  %s %s
EOF
''' % (opt_file, r"%chk", chk_file, keyword, molname + td_suffix, charge, multiplicity))

    os.system('''awk '{if (NR>5 && $1 !~ "[0-9]") print }' tmp.gjf >> %s''' % opt_file)
    os.remove('tmp.gjf')

    return opt_file, chk_file


def generate_energy(original_st: str, charge: int, multiplicity: int, dft: str = 'B3LYP', basis_set: str = '6-31g*', solvent: str = 'water', correct: bool = True, td: bool = False):
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
    if os.path.exists(energy_file):
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

    os.system('''awk '{if (NR>5 && $1 !~ "[0-9]") print }' tmp.gjf >> %s''' % energy_file)
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
        logger.warning(
            'You are using gaussian 09, that may cause some unknown errors.\nGaussian 16 is recommend.')
        gaussian = g09
    else:
        logger.error('Gaussian is not installed.')
        return None
    return gaussian


def system_default(gauss_path: str, cpu_count: int, memory: str):
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
    with open(default_route, 'w') as f:
        f.write('-P- %s\n' % cpu_count)
        f.write('-M- %s\n' % memory)
    logger.debug('Default system setting changed: CPU = %s  Mem = %s' %
                 (cpu_count, memory))


def generate_fchk(chk_file: str):
    '''
    生成fchk文件

    Return
    ----------
    str
        生成的fchk文件
    '''
    fchk_file = chk_file.split('.')[0] + '.fchk'
    os.system('formchk %s > /dev/null' % chk_file)
    logger.debug('fchk file %s is saved.' % fchk_file)
    return fchk_file

def _check_gauss_finished(line: str):
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

def _get_system_info(gauss_path: str):
    '''
    读取当前高斯计算资源文件设定

    Parameter
    ---------
    gauss_path : str
        高斯可执行文件路径
    Return
    ---------
    tuple[str, str]
        CPU, memory
    '''

    default_route = os.path.dirname(gauss_path) + '/Default.Route'
    if not os.path.exists(default_route):
        with open(default_route, 'w') as f:
            f.write('-P- 8\n')
            f.write('-M- 4GB')
    with open(default_route, 'r') as f:
        cpu_info, mem_info = f.read().splitlines()

    return cpu_info[4:], mem_info[4:]


def get_mo(fchk_file: str):
    '''
    获取HOMO/LUMO分子轨道信息

    Parameters
    ----------
    fchk_file : str
        高斯计算检查点文件(非二进制)

    Return
    ----------
    dict[str, dict]
        HOMO, LUMO分子轨道相关信息, gap值  

        {
        'homo': 
                {'index': homo_index, 
                'energy': homo_energy}, 
                
        'lumo': 
            {'index': lumo_index, 
            'energy': lumo_energy}, 

        'gap': gap value
        }
    '''
    mo_pipe = subprocess.Popen('Multiwfn %s ' % fchk_file, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    mo_info = ''

    for line in iter(mo_pipe.stdout.readline, b''):
        if re.search(r'Ghost atoms', line.decode('utf-8')):
            logger.warning('Ghost atoms (Bq) are found in this file')
            mo_pipe.stdin.write(b'y\n')
            mo_pipe.stdin.flush()
        if re.search(r'300 Other functions', line.decode('utf-8')):
            mo_pipe.stdin.write(b'0\n')
            mo_pipe.stdin.write(b'q\n')
            mo_pipe.stdin.flush()
        mo_info += line.decode('utf-8')
    mo_pipe.stdin.close()
    mo_pipe.stdout.close()
        
    homo_index = int(re.search(r'[ \d]+(?=is HOMO)', mo_info).group().strip())
    homo_energy = '%.6f' % (float(re.search(r'(?<=is HOMO, energy:)[-\d. ]+', mo_info).group().strip()) * RATIO)
    lumo_index = int(re.search(r'[ \d]+(?=is LUMO)', mo_info).group().strip())
    lumo_energy = '%.6f' % (float(re.search(r'(?<=is LUMO, energy:)[-\d. ]+', mo_info).group().strip()) * RATIO)
    gap = '%.6f' % (float(re.search(r'(?<=HOMO-LUMO gap:)[-\d. ]+', mo_info).group().strip()) * RATIO)

    return {'homo': {'index': homo_index, 'energy': homo_energy}, 'lumo': {'index': lumo_index, 'energy': lumo_energy}, 'gap': gap}


def cube_file_generate(fchk_file: str, mo: int):
    '''
    生成分子轨道cube Grid文件

    Parameters
    ----------
    fchk_file : str
        高斯计算检查点文件(非二进制)
    mo : int
        分子轨道(MO)编号

    Return
    ----------
    str
        生成的Grid文件名
    '''

    logger.debug('Extracting MO %s from %s' % (mo, fchk_file))
    _excute_pipe = subprocess.Popen('Multiwfn %s ' % fchk_file, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

    for line in iter(_excute_pipe.stdout.readline, b''):
        if re.search(r'Ghost atoms', line.decode('utf-8')):
            _excute_pipe.stdin.write(b'y\n')
            _excute_pipe.stdin.flush()
        if re.search(r'300 Other functions', line.decode('utf-8')):
            _excute_pipe.stdin.write(b'200\n')
            _excute_pipe.stdin.write(b'3\n')
            _excute_pipe.stdin.write(b'%s\n' % bytes(str(mo), 'utf-8'))
            _excute_pipe.stdin.write(b'2\n')
            _excute_pipe.stdin.write(b'1\n')
            _excute_pipe.stdin.write(b'0\n')
            _excute_pipe.stdin.write(b'q\n')
            _excute_pipe.stdin.flush()

    return 'orb' + str(mo).rjust(6, '0') + '.cub'

class Daemon:
    def __init__(self, cmd, pidfile='/tmp/daemon.pid', stdin='/dev/null', stdout='/dev/null', stderr='/dev/null'):
        self.stdin = stdin
        self.stdout = stdout
        self.stderr = stderr
        self.pidfile = pidfile
        self.cmd = cmd

    def daemonize(self):
        if os.path.exists(self.pidfile):
            raise RuntimeError('Already running.')
        pid = os.fork()
        # First fork (detaches from parent)
        try:
            if pid > 0:
                raise SystemExit(0)
        except OSError as e:
            raise RuntimeError('fork #1 faild: {0} ({1})\n'.format(e.errno, e.strerror))

        #os.chdir('/')
        os.setsid()
        os.umask(0o22)

        # Second fork (relinquish session leadership)
        _pid = os.fork()
        try:
            if _pid > 0:
                raise SystemExit(0)
        except OSError as e:
            raise RuntimeError('fork #2 faild: {0} ({1})\n'.format(e.errno, e.strerror))

        # Flush I/O buffers
        sys.stdout.flush()
        sys.stderr.flush()

        # Replace file descriptors for stdin, stdout, and stderr
        '''        
        with open(self.stdin, 'rb', 0) as f:
            os.dup2(f.fileno(), sys.stdin.fileno())
        with open(self.stdout, 'ab', 0) as f:
            os.dup2(f.fileno(), sys.stdout.fileno())
        with open(self.stderr, 'ab', 0) as f:
            os.dup2(f.fileno(), sys.stderr.fileno())
        '''
        # Write the PID file
        with open(self.pidfile, 'w') as f:
            print(os.getpid(), file=f)
        
        # Arrange to have the PID file removed on exit/signal
        atexit.register(lambda: os.remove(self.pidfile))

        signal.signal(signal.SIGTERM, self.__sigterm_handler)
        

    # Signal handler for termination (required)
    @staticmethod
    def __sigterm_handler(signo, frame):
        raise SystemExit(1)

    def start(self):
        try:
            self.daemonize()
        except RuntimeError as e:
            print(e, file=sys.stderr)
            raise SystemExit(1)

        self.run()

    def stop(self):
        try:
            if os.path.exists(self.pidfile):
                with open(self.pidfile) as f:
                    os.kill(int(f.read()), signal.SIGTERM)
            else:
                print('Not running.', file=sys.stderr)
                raise SystemExit(1)
        except OSError as e:
            if 'No such process' in str(e) and os.path.exists(self.pidfile): 
                os.remove(self.pidfile)

    def restart(self):
        self.stop()
        self.start()

    def run(self):
        subprocess.Popen(self.cmd, shell=True, bufsize=1)
