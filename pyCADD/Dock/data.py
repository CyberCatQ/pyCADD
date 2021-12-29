import csv
import logging
import os
import re

from pyCADD.Dock.prepare import load_st
from pyCADD.utils.check import check_file
from schrodinger import structure as struc
logger = logging.getLogger('pyCADD.data')

def extra_data(file_path:str) -> dict:
    '''
    从对接或计算完成的Maestro文件中提取一般数据

    Parameters
    ----------
    file_path : str
        对接完成的文件PATH

    Return
    ----------
    dict
        Properties : Values

        '''
    logger.debug('Prepare to extract data from file %s' % file_path)
    file = os.path.basename(file_path).split('.')[0]
    if not check_file(file_path):
        return

    # 从文件名获取信息
    pdbid = file.split('_')[0]
    original_lig = file.split('_')[1]
    ligname = file.split('_')[4]
    precision = file.split('_')[5]
    
    st = struc.StructureReader(file_path)
    pro_st = next(st)
    lig_st = next(st)

    prop_dic = {}

    # 需要提取的Property 公共项
    prop_dic['PDB'] = pdbid
    prop_dic['Ligand'] = ligname
    prop_dic['Original'] = original_lig
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

    sitemap_file = '%s_sitemap_out.maegz' % pdbid

    if os.path.exists(sitemap_file):
        site_st = load_st(sitemap_file)
        try:
            prop_dic['Site_Score'] = site_st.property['r_sitemap_SiteScore']
            prop_dic['Volume'] = site_st.property['r_sitemap_volume']
        except KeyError:
            pass
    
    return prop_dic

def save_data(data_dic:dict, pdbid:str, ligname:str, precision:str='SP', additional_col:list=[]):
    '''
    储存一般数据为csv

    Parameter
    ----------
    data_dic : dict
        数据内容 {property : data}
    pdbid : str
        PDB ID
    ligname : str
        配体ID(RCSB ID)
    precision : str
        对接精度
    additional_col : list
        额外需要提取的数据列名(maestro原始名称)
    '''
    prop_xp = ['PDB', 'Ligand', 'Original', 'Docking_Score', 'MMGBSA_dG_Bind', 'rmsd', 'precision', 'Site_Score', 'Volume','ligand_efficiency', 
            'XP_Hbond', 'rotatable_bonds', 'ecoul', 'evdw', 'emodel', 'energy', 'einternal', 'activity']
    prop_sp = ['PDB', 'Ligand', 'Original', 'Docking_Score', 'MMGBSA_dG_Bind', 'rmsd', 'precision', 'Site_Score', 'Volume', 'ligand_efficiency', 
            'rotatable_bonds', 'ecoul', 'evdw', 'emodel', 'energy', 'einternal','lipo', 'hbond', 'metal', 'rewards', 'erotb', 'esite', 'activity']
    
    if precision == 'SP':
        fields = prop_sp
    elif precision == 'XP':
        fields = prop_xp
    
    if additional_col:
        for col in additional_col:
            fields.append(col)
        
    with open(pdbid + '_FINAL_RESULTS_%s.csv' % ligname, 'w', encoding='UTF-8', newline='') as f:

        writer = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()            # 写入标头
        writer.writerows([data_dic])    # 写入数据 自动匹配标头列
    
    logger.debug('%s-%s Data Saved.\n' % (pdbid, ligname))

def extra_admet_data(admet_file_path:str) -> dict:

    '''
    提取ADMET数据

    Parameter
    ----------
    admet_file_path : str
        ADMET计算结果文件PATH

    Return
    ----------
    dict
        提取数据组成的字典
    '''
    logger.debug('\nPrepare to extract ADMET data from file %s' % admet_file_path)

    admet_dic = {}
    st = load_st(admet_file_path)
    admet_file = os.path.basename(admet_file_path)
    ligname = admet_file.split('.')[0].split('_')[-2]

    admet_dic['Ligand'] = ligname

    for _key in st.property.keys():
        # Qikprop生成的分子描述符键名
        if re.match(r'^[ri]_qp_.+', _key):                           
            key = re.search(r'(?<=[ri]_qp_).+', _key).group()
            admet_dic[key] = st.property[_key]
        
    return admet_dic

def save_admet_data(admet_dic:dict, ligname:str):
    '''
    保存ADMET数据文件

    Parameters
    ----------
    admet_dic : dic
        ADMET数据dict
    ligname : str
        配体ID(RCSB ID)
    '''

    prop_admet = list(admet_dic.keys())

    with open('%s_ADMET.csv' % ligname, 'w', encoding='UTF-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=prop_admet)
        writer.writeheader()                # 写入标头
        writer.writerows([admet_dic])       # 写入数据 自动匹配标头列

    logger.debug('ADMET data of %s saved.\n' % ligname)