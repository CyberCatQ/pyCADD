import csv
import logging
import re
import os
import pandas as pd
from typing import List

from pyCADD.Dock.common import DockResultFile, LigandFile
from pyCADD.Dock.config import DefaultDataConfig, DataConfig
logger = logging.getLogger(__name__)

def _save_data(output_file:str, data:list, fields:list):
    '''
    储存数据为csv

    Parameter
    ----------
    output_file : str
        输出文件路径
    data_dic : dict
        数据内容 {property : data}
    fields : List[str]
        数据提取项配置
    '''
    with open(output_file, 'w', encoding='UTF-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()            # 写入标头
        writer.writerows(data)          # 写入数据 自动匹配标头列

def extra_docking_data(dock_result_file:DockResultFile) -> dict:
    '''
    从对接或计算完成的Maestro文件中提取一般数据或配体文件

    Parameters
    ----------
    dock_result_file : DockResultFile
        对接完成的文件

    Return
    ----------
    dict
        数据Properties : Values
    '''
    logger.debug(f'Prepare to extract data from file {dock_result_file.file_path}')
    pdbid = dock_result_file.pdbid
    internal_ligand_name = dock_result_file.internal_ligand_name
    docking_ligand_name = dock_result_file.docking_ligand_name
    precision = dock_result_file.precision
    
    lig_st = dock_result_file.docking_ligand_st

    prop_dic = {}

    # 需要提取的Property 公共项
    prop_dic['PDB'] = pdbid
    prop_dic['Ligand'] = docking_ligand_name
    prop_dic['Original'] = internal_ligand_name
    prop_dic['Docking_Score'] = lig_st.property['r_i_docking_score']  # 对接分数
    prop_dic['rotatable_bonds'] = lig_st.property['i_i_glide_rotatable_bonds']
    prop_dic['ligand_efficiency'] = lig_st.property['r_i_glide_ligand_efficiency']
    prop_dic['evdw'] = lig_st.property['r_i_glide_evdw']
    prop_dic['ecoul'] = lig_st.property['r_i_glide_ecoul']
    prop_dic['energy'] = lig_st.property['r_i_glide_energy']
    prop_dic['einternal'] = lig_st.property['r_i_glide_einternal']
    prop_dic['emodel'] = lig_st.property['r_i_glide_emodel']
    prop_dic['precision'] = precision

    # 其他属性
    for key in lig_st.property.keys():
        prop_dic[key] = lig_st.property[key]

    # rmsd项
    try:
        prop_dic['rmsd'] = lig_st.property['r_i_glide_rmsd_to_input']
    except KeyError:
        pass
    # MMGBSA结合能项
    try:
        prop_dic['MMGBSA_dG_Bind'] = lig_st.property['r_psp_MMGBSA_dG_Bind']
    except KeyError:
        pass
    # 活性标签
    try:
        prop_dic['activity'] = lig_st.property['b_user_Activity']
    except KeyError:
        pass
    # SP精度特有项
    if precision == 'SP':
        prop_dic['lipo'] = lig_st.property['r_i_glide_lipo']
        prop_dic['hbond'] = lig_st.property['r_i_glide_hbond']
        prop_dic['metal'] = lig_st.property['r_i_glide_metal']
        prop_dic['rewards'] = lig_st.property['r_i_glide_rewards']
        prop_dic['erotb'] = lig_st.property['r_i_glide_erotb']
        prop_dic['esite'] = lig_st.property['r_i_glide_esite']
    #XP精度特有项
    elif precision == 'XP':
        prop_dic['XP_Hbond'] = lig_st.property['r_glide_XP_HBond']

    return prop_dic

def save_docking_data(dock_result_file:DockResultFile, configs:DataConfig=None):
    '''
    储存一般数据为csv

    Parameter
    ----------
    data_dic : dict
        数据内容 {property : data}
    configs : DataConfig
        数据提取项配置
    '''
    data_dic = extra_docking_data(dock_result_file)
    pdbid = data_dic['PDB']
    ligand_name = data_dic['Ligand']
    precision = data_dic['precision']

    property_config = DefaultDataConfig(precision) if configs is None else configs
    fields = property_config.properties

    output_file = f'{pdbid}_{ligand_name}_FINAL_RESULTS.csv'
    _save_data(output_file, [data_dic], fields)
    logger.debug('%s-%s Data Saved.' % (pdbid, ligand_name))
    return output_file

def save_redocking_data(data_list:list, precision:str='SP', configs:DataConfig=None, save_dir:str=None) -> None:
    '''
    储存回顾性Ensemble Docking对接数据

    Parameter
    ----------
    data_list : list
        回顾性Ensemble Docking对接数据
    precision : str
        对接精度
    configs : DataConfig
        数据提取项配置
    save_dir : str
        储存目录
    '''
    save_dir = os.getcwd() if save_dir is None else save_dir
    property_config = DefaultDataConfig(precision) if configs is None else configs
    fields = property_config.properties
    reference_results_file_path = os.path.join(save_dir, 'ENSEMBLE_DOCKING_RESULTS_REF.csv')

    redock_df = pd.DataFrame(data_list)
    # redock的内源配体默认为阳性
    redock_df['activity'] = 1
    redock_df.to_csv(reference_results_file_path, index=False)
    # Docking Score矩阵生成
    ligand_df_group = redock_df.groupby('Ligand')
    total_redock_data = []
    total_rmsd_data = []

    for ligand_name, ligand_df in ligand_df_group:
        total_redock_data.append({'Ligand': ligand_name, **ligand_df[['PDB', 'Docking_Score']].set_index('PDB').loc[:,'Docking_Score'].to_dict()})
        total_rmsd_data.append({'Ligand': ligand_name, **ligand_df[['PDB', 'rmsd']].set_index('PDB').loc[:,'rmsd'].to_dict()})

    docking_score_matrix_df = pd.DataFrame(total_redock_data)
    rmsd_matrix_df = pd.DataFrame(total_rmsd_data)
    docking_score_matrix_df.to_csv(os.path.join(save_dir, 'reference_matrix.csv'), index=False)
    rmsd_matrix_df.to_csv(os.path.join(save_dir, 'reference_rmsd_matrix.csv'), index=False)
    
def save_ensemble_docking_data(data_list:list, precision:str='SP', configs:DataConfig=None, save_dir:str=None) ->None:
    '''
    分类储存Ensemble Docking对接数据

    Parameter
    ----------
    data_list : list
        Ensemble Docking对接数据列表
    precision : str
        对接精度
    configs : DataConfig
        数据提取项配置
    '''
    save_dir = os.getcwd() if save_dir is None else save_dir
    property_config = DefaultDataConfig(precision) if configs is None else configs
    fields = property_config.properties
    final_results_file_path = os.path.join(save_dir, 'ENSEMBLE_DOCKING_RESULTS.csv')
    total_df = pd.DataFrame(data_list)
    total_df.to_csv(final_results_file_path, index=False)

    # 分类储存
    total_df_group = total_df.groupby('PDB')
    for pdbid, pdb_df in total_df_group:
        pdb_df = pd.DataFrame(pdb_df, columns=fields)
        pdb_df.to_csv(os.path.join(save_dir, f'{pdbid}_ENSEMBLE_DOCKING_RESULTS.csv'), index=False)
    
    # Docking Score矩阵生成
    ligand_df_group = total_df.groupby('Ligand')
    total_ligand_data = []
    for ligand_name, ligand_df in ligand_df_group:
        total_ligand_data.append({'Ligand': ligand_name, **ligand_df[['PDB', 'Docking_Score']].set_index('PDB').loc[:,'Docking_Score'].to_dict()})
    matrix_df = pd.DataFrame(total_ligand_data)
    matrix_df.to_csv(os.path.join(save_dir, 'matrix.csv'), index=False)
        
def extra_admet_data(admet_file:LigandFile) -> List[dict]:

    '''
    提取ADMET数据

    Parameter
    ----------
    admet_file : LigandFile
        要提取的ADMET计算结果文件

    Return
    ----------
    List[dict]
        数据Properties : Values
    '''
    logger.debug(f'Prepare to extract ADMET data from file {admet_file.file_path}')

    admet_list = []

    for index, st in enumerate(admet_file.st_reader):
        curr_prop_dict = {}
        curr_prop_dict['Title'] = st.property['s_m_title']
        curr_prop_dict['index'] = index
        for _key in st.property.keys():
            # Qikprop生成的分子描述符键名
            _match = re.match(r'^[ri]_qp_.+', _key)
            if _match:
                curr_prop_dict[_match.group()] = st.property[_key]
        admet_list.append(curr_prop_dict)
        
    return admet_list

def save_admet_data(admet_file:LigandFile):
    '''
    保存ADMET数据文件

    Parameters
    ----------
    admet_file : LigandFile
        要保存的ADMET计算结果文件
    '''
    data_list = extra_admet_data(admet_file)
    admet_property = list(data_list[0].keys())
    output_file = f'{admet_file.file_prefix}_ADMET_RESULT.csv'

    _save_data(output_file, data_list, admet_property)
    logger.debug(f'ADMET data {output_file} saved.')
    return output_file