from pyCADD.Dock.core import launch
from pyCADD.utils.check import check_file
from pyCADD.utils.getinfo import get_config, get_gridfile_path_list, get_project_dir

import os
import logging

logger = logging.getLogger('pyCADD.VSW.core')

def read_gene_config(gene_config_path):
    '''
    从文件中读取可供筛选的受体基因

    Parameter
    ----------
    gene_config_path : str
        受体信息配置文件路径

    Return
    ----------
    dict
        受体配置信息
    '''

    return dict(get_config(gene_config_path)['GENE'])

def read_database_config(database_config_path):
    '''
    从文件中读取化合物库路径信息

    Parameter
    ----------
    database_config_path : str
        化合物库路径信息文件路径
    
    Return
    ----------
    dict
        化合物库路径配置信息
    '''

    return dict(get_config(database_config_path)['DATABASE'])

def gen_input_file(recep_list:list, lig_file:str, jobname:str=''):
    '''
    生成vsw输入文件

    Parameter
    ---------
    recep_list : list
        将要作为受体进行VSW的所有受体(PDBID, 配体ID)组成的列表
    lig_file : str
        用于VSW的化合物库文件PATH

    Return
    ----------
    str
        生成的inp输入文件名
    '''
    cwd = get_project_dir()
    vsw_dir = cwd + '/vsw/'
    
    grid_path_list = get_gridfile_path_list(recep_list)                         # 生成Grid文件名列表
    os.chdir(vsw_dir)
    grid_num = len(grid_path_list)                                              # 即将为输入文件写入的Pipeline数量 即PDB总数

    input_file = jobname + '.inp'
    with open(input_file, 'w') as f:                                            # VSW 输入文件编写
        f.write('# Run as: $SCHRODINGER/vsw <inputfile> \n')

        # 载入配体数据集
        f.write('\n[SET:ORIGINAL_LIGANDS]\n')                                     
        f.write('    VARCLASS   Structures\n')
        f.write('    FILES   %s,\n' % lig_file)                                 # 化合物库PATH

        # 载入所有受体数据集
        for grid_path in grid_path_list:
            index = grid_path_list.index(grid_path) + 1                         # 默认从0开始索引
            f.write('\n[SET:GRID_%s]\n' % index)
            f.write('    VARCLASS   Grid\n')
            f.write('    FILE   %s\n' % grid_path)
            
        # 执行步骤
        for n in range(1, grid_num + 1):

            # PRE_DOCK_HTVS 参数                                   
            f.write('\n[STAGE:PRE_DOCK_HTVS_%s]\n' % n)                         
            f.write('    STAGECLASS   gencodes.RecombineStage\n')
            f.write('    INPUTS   ORIGINAL_LIGANDS,\n')
            f.write('    OUTPUTS   RECOMBINE_OUT_%s,\n' % n)
            f.write('    NUMOUT   njobs\n')
            f.write('    OUTFORMAT   maegz\n')
            f.write('    MIN_SUBJOB_STS   4000\n')
            f.write('    MAX_SUBJOB_STS   40000\n')
            f.write('    GENCODES   YES\n')
            f.write('    OUTCOMPOUNDFIELD   s_vsw_compound_code\n')
            f.write('    OUTVARIANTFIELD   s_vsw_variant\n')
            f.write('    UNIQUEFIELD   NONE\n')

            # DOCK_HTVS 参数
            f.write('\n[STAGE:DOCK_HTVS_%s]\n' % n)
            f.write('    STAGECLASS   glide.DockingStage\n')
            f.write('    INPUTS   RECOMBINE_OUT_%s, GRID_%s\n' % (n, n))
            f.write('    OUTPUTS   HTVS_OUT_%s,\n' % n)
            f.write('    RECOMBINE   NO\n')
            f.write('    PRECISION   HTVS\n')
            f.write('    UNIQUEFIELD   s_vsw_compound_code\n')
            f.write('    PERCENT_TO_KEEP   10.0\n')                 # 保留前10%
            f.write('    DOCKING_METHOD   confgen\n')
            f.write('    POSES_PER_LIG   1\n')
            f.write('    BEST_BY_TITLE   NO\n')
            f.write('    LIG_VSCALE   0.8\n')
            f.write('    LIG_CCUT   0.15\n')
            f.write('    MAXATOMS   300\n')
            f.write('    MAXROTBONDS   50\n')
            f.write('    AMIDE_MODE   penal\n')
            f.write('    POSE_OUTTYPE   LIB\n')
            f.write('    POSTDOCK   YES\n')
            f.write('    POSTDOCKSTRAIN   NO\n')
            f.write('    COMPRESS_POSES   YES\n')
            f.write('    EPIK_PENALTIES   YES\n')
            f.write('    FORCEPLANAR   NO\n')

            # PULL_HTVS 参数
            f.write('\n[STAGE:PULL_HTVS_%s]\n' % n)
            f.write('    STAGECLASS   pull.PullStage\n')
            f.write('    INPUTS   HTVS_OUT_%s, RECOMBINE_OUT_%s\n' % (n, n))
            f.write('    OUTPUTS   HTVS_OUT_ORIG_%s,\n' % n)
            f.write('    UNIQUEFIELD   s_vsw_compound_code\n')

            # PRE_DOCK_SP 参数
            f.write('\n[STAGE:PRE_DOCK_SP_%s]\n' % n)
            f.write('    STAGECLASS   gencodes.RecombineStage\n')
            f.write('    INPUTS   HTVS_OUT_ORIG_%s,\n' % n)
            f.write('    OUTPUTS   DOCK_SP_%s_INPUT,\n' % n)
            f.write('    NUMOUT   njobs\n')
            f.write('    OUTFORMAT   maegz\n')
            f.write('    MIN_SUBJOB_STS   300\n')
            f.write('    MAX_SUBJOB_STS   3000\n')
            f.write('    GENCODES   NO\n')
            f.write('    UNIQUEFIELD   s_vsw_compound_code\n')

            # DOCK_SP 参数
            f.write('\n[STAGE:DOCK_SP_%s]\n' % n)
            f.write('    STAGECLASS   glide.DockingStage\n')
            f.write('    INPUTS   DOCK_SP_%s_INPUT, GRID_%s\n' % (n, n))
            f.write('    OUTPUTS   SP_OUT_%s,\n' % n)
            f.write('    RECOMBINE   NO\n')
            f.write('    PRECISION   SP\n')
            f.write('    UNIQUEFIELD   s_vsw_compound_code\n')
            f.write('    PERCENT_TO_KEEP   50\n')                              # 保留前50%
            f.write('    DOCKING_METHOD   confgen\n')
            f.write('    POSES_PER_LIG   1\n')
            f.write('    WRITE_XP_DESC   NO\n')
            f.write('    NENHANCED_SAMPLING   1\n')
            f.write('    BEST_BY_TITLE   NO\n')
            f.write('    LIG_VSCALE   0.8\n')
            f.write('    LIG_CCUT   0.15\n')
            f.write('    MAXATOMS   300\n')
            f.write('    MAXROTBONDS   50\n')
            f.write('    AMIDE_MODE   penal\n')
            f.write('    POSE_OUTTYPE   LIB\n')
            f.write('    POSTDOCK   YES\n')
            f.write('    POSTDOCKSTRAIN   NO\n')
            f.write('    COMPRESS_POSES   YES\n')
            f.write('    EPIK_PENALTIES   YES\n')
            f.write('    FORCEPLANAR   NO\n')

            # PULL_SP 参数
            f.write('\n[STAGE:PULL_SP_%s]\n' % n)
            f.write('    STAGECLASS   pull.PullStage\n')
            f.write('    INPUTS   SP_OUT_%s, HTVS_OUT_ORIG_%s\n' % (n, n))
            f.write('    OUTPUTS   SP_OUT_ORIG_%s,\n' % n)
            f.write('    UNIQUEFIELD   s_vsw_variant\n')

            # PRE_DOCK_XP 参数
            f.write('\n[STAGE:PRE_DOCK_XP_%s]\n' % n)
            f.write('    STAGECLASS   gencodes.RecombineStage\n')
            f.write('    INPUTS   SP_OUT_ORIG_%s,\n' % n)
            f.write('    OUTPUTS   DOCK_XP_%s_INPUT,\n' % n)
            f.write('    NUMOUT   njobs\n')
            f.write('    OUTFORMAT   maegz\n')
            f.write('    MIN_SUBJOB_STS   20\n')
            f.write('    MAX_SUBJOB_STS   200\n')
            f.write('    GENCODES   NO\n')
            f.write('    UNIQUEFIELD   s_vsw_compound_code\n')

            # DOCK_XP 参数
            f.write('\n[STAGE:DOCK_XP_%s]\n' % n)
            f.write('    STAGECLASS   glide.DockingStage\n')
            f.write('    INPUTS   DOCK_XP_%s_INPUT, GRID_%s\n' % (n, n))
            f.write('    OUTPUTS   XP_OUT_%s,\n' % n)
            f.write('    RECOMBINE   NO\n')
            f.write('    PRECISION   XP\n')
            f.write('    UNIQUEFIELD   s_vsw_compound_code\n')
            f.write('    PERCENT_TO_KEEP   20.0\n')                           # 保留前20%
            f.write('    DOCKING_METHOD   confgen\n')
            f.write('    POSES_PER_LIG   1\n')
            f.write('    WRITE_XP_DESC   NO\n')
            f.write('    BEST_BY_TITLE   YES\n')
            f.write('    LIG_VSCALE   0.8\n')
            f.write('    LIG_CCUT   0.15\n')
            f.write('    MAXATOMS   300\n')
            f.write('    MAXROTBONDS   50\n')
            f.write('    AMIDE_MODE   penal\n')
            f.write('    POSE_OUTTYPE   PV\n') 
            f.write('    POSTDOCK   YES\n')
            f.write('    POSTDOCKSTRAIN   NO\n')
            f.write('    COMPRESS_POSES   YES\n')
            f.write('    EPIK_PENALTIES   YES\n')
            f.write('    FORCEPLANAR   NO\n')

            # MMGBSA 参数
            f.write('\n[STAGE:MMGBSA_%s]\n' % n)
            f.write('    STAGECLASS   prime.MMGBSAStage\n')
            f.write('    INPUTS   XP_OUT_%s,\n' % n)
            f.write('    OUTPUTS   MMGBSA_%s,\n' % n)
        
        # 数据合并
        f.write('\n[STAGE:DOCKMERGE]\n')
        f.write('    STAGECLASS   glide.MergeStage\n')
        f.write('    INPUTS   ')
        for n in range(1, grid_num + 1):
            f.write('MMGBSA_%s,' % n)
        f.write('\n')
        f.write('    OUTPUTS   FINAL_DOCK_OUT,\n')
        f.write('    MAXPERLIG   1\n')
        f.write('    NREPORT   10000\n')
        f.write('    UNIQUEFIELD   s_vsw_compound_code\n')
        f.write('    OFFSETS   ')
        for n in range(1, grid_num + 1):
            f.write('0.0, ')
        f.write('\n')
        f.write('    MARKERS   ')
        for n in range(1, grid_num + 1):
            f.write('%s, ' % n)
        f.write('\n')

        # 数据输出
        f.write('\n[USEROUTS]\n')
        f.write('    USEROUTS   ')
        for n in range(1, grid_num + 1):
            f.write('XP_OUT_%s, MMGBSA_%s, ' % (n, n))
        f.write('FINAL_DOCK_OUT\n')
        f.write('    STRUCTOUT   FINAL_DOCK_OUT')

    os.chdir(cwd)
    return input_file