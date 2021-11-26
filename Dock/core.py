import os

from pyCADD.Dock.prepare import load_st
from schrodinger.job import jobcontrol as jc


def launch(cmd):
    '''
        使用jobcontrol启动一项job并等待结束

        Parameters
        ----------
        cmd ： str
            等待执行的命令字符串
        '''

    cmd_list = cmd.split(' ')  # launch_job以列表形式提交参数
    job = jc.launch_job(cmd_list)
    print('JobId: %s' % job.JobId, end='\n')
    print('Job Name: %s' % job.Name, end='\n')
    print('Job Status: %s' % job.Status, end='\n')
    job.wait()  # 阻塞进程 等待Job结束

def keep_chain(pdbfile:str, chain_name:str) -> str:
    '''
    读取PDB晶体文件并将单一链的结构输出为pdb文件

    Parameters
    ----------
    pdbfile : str
        待处理原始文件PATH
    chain_name : str
        要保留的链名称

    Return
    ----------
    str
        保留单链结构的文件名
    '''

    st = load_st(pdbfile)  # 读取原始PDB结构
    st_chain_only = st.chain[chain_name].extractStructure()
    singlechain_file = '%s_chain_%s.mae' % (pdbfile.split('.')[0], chain_name)
    st_chain_only.write(singlechain_file)
    return singlechain_file

def minimize(pdbfile:str) -> str:
    '''
    调用prepwizard模块自动优化PDB结构文件

    Parameters
    ----------
    pdbfile : str
        需要优化的文件PATH
    
    Return
    ----------
    str
        完成优化后的文件名

    '''

    minimized_file = str(pdbfile.split('.')[0] + '_minimized.mae')

    if os.path.exists(minimized_file):  # 如果已经进行过优化 为提高效率而跳过优化步骤
        return minimized_file

    prepwizard_command = 'prepwizard -f 3 -r 0.3 -propka_pH 7.0 -disulfides -s -j %s-Minimize %s %s' % (pdbfile.split('.')[0],
                                                                                                            pdbfile, minimized_file)
    launch(prepwizard_command)   # 阻塞至任务结束

    # 判断Minimized任务是否完成(是否生成Minimized结束的结构文件)
    if not os.path.exists(minimized_file):  
        # 无法被优化的晶体结构
        raise RuntimeError('%s Crystal Minimization Process Failed' % pdbfile.split('.')[0])  
    else:
        print('\nPDB Minimized File', minimized_file, 'Saved.\n')
        return minimized_file
