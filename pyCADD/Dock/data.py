import csv
import logging
import re
import os
import pandas as pd
from typing import List
import xlsxwriter
from io import BytesIO
from rdkit.Chem import rdDepictor, PandasTools, MolToSmiles, Draw, MolFromSmiles
rdDepictor.SetPreferCoordGen(True)

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
    # fields = property_config.properties
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

# Quick Report Functions
def _get_smiles(st_path):
    '''
    获取SMILES
    '''
    from schrodinger.structure import StructureReader
    from schrodinger import adapter

    if not os.path.exists(st_path):
        logger.warning(f'{st_path} not exists.')
        return None
    _st = StructureReader.read(st_path)
    _mol = adapter.to_rdkit(_st)
    _smiles = MolToSmiles(_mol)
    return _smiles

def _generate_img(dataframe, ligand_col:str, ligand_save_dir:str):
    '''
    为参考数据生成结构图

    Parameter
    ----------
    dataframe : pd.DataFrame
        参考数据
    ligand_col : str
        参考数据中的Ligand列名
    ligand_save_dir : str
        参考结构结构保存路径
    '''
    
    _smiles_list = []
    for ligand_name in dataframe[ligand_col]:
        st_path = os.path.join(ligand_save_dir, f'{ligand_name}.mae')
        _smiles = _get_smiles(st_path)
        _smiles_list.append({ligand_col: ligand_name, 'SMILES': _smiles})
    smiles_df = pd.DataFrame(_smiles_list)
    PandasTools.AddMoleculeColumnToFrame(smiles_df, 'SMILES', 'Structure')
    smiles_df.drop('SMILES', axis=1, inplace=True)
    dataframe = dataframe.merge(smiles_df, on=ligand_col, how='left')
    dataframe = dataframe[['Receptor', 'PDB', ligand_col, 'Structure', 'Docking_Score', 'rmsd']]
    return dataframe

def _write_xlsx(df, file_path, info_dict:dict=None):
    '''
    将含有结构图的dataframe写入文件

    Parameters
    ----------
    df : pd.DataFrame
        含有结构图的dataframe
    file_path : str
        要写入的文件路径
    '''
    info_dict = {} if info_dict is None else info_dict

    output_workbook = xlsxwriter.Workbook(file_path)
    output_worksheet = output_workbook.add_worksheet()
    
    base_fmt_dict = {'font_name': u'等线', 'font_color': 'black', 'align': 'center', 'valign': 'vcenter'}
    title_fmt = output_workbook.add_format({**base_fmt_dict, 'bold': True, 'font_size': 14, 'bottom':1})
    subtitle_fmt = output_workbook.add_format({**base_fmt_dict, 'bold': True, 'font_size': 12, 'top': 1, 'bottom': 1})
    center_fmt = output_workbook.add_format({**base_fmt_dict, 'num_format': '0.0000', 'font_size': 12})
    column_fmt = output_workbook.add_format({**base_fmt_dict, 'font_size': 12})
    bold_fmt = output_workbook.add_format({**base_fmt_dict, 'bold': True, 'font_size': 12})
    highlight_fmt = output_workbook.add_format({**base_fmt_dict, 'font_size': 12, 'num_format': '0.0000', 'fg_color': 'green'})
    endline_fmt = output_workbook.add_format({**base_fmt_dict, 'font_size': 12, 'top': 1})

    ligand_name = info_dict.get('ligand_name')
    title = f'{ligand_name} Docking Score Report'
    output_worksheet.merge_range('A1:G1', title, title_fmt)

    output_worksheet.set_column('A:B', width=10, cell_format=column_fmt)
    output_worksheet.set_column('C:C', width=16, cell_format=column_fmt)
    output_worksheet.set_column('D:D', width=20.91, cell_format=column_fmt)
    output_worksheet.set_column('E:E', width=14, cell_format=column_fmt)
    output_worksheet.set_column('F:F', width=10, cell_format=column_fmt)
    output_worksheet.set_column('G:G', width=14, cell_format=column_fmt)

    output_worksheet.set_row(1, 80)
    output_worksheet.set_row(2, 80)
    output_worksheet.set_row(3, 20)

    output_worksheet.write('A2', info_dict.get('higher_count'))
    output_worksheet.write('B2', 'Higher', bold_fmt)
    output_worksheet.write('A3', info_dict.get('lower_count'))
    output_worksheet.write('B3', 'Lower', bold_fmt)

    output_worksheet.merge_range('C2:C3', info_dict.get('ligand_name'), center_fmt)
    _imgdata = BytesIO()
    _img = Draw.MolToImage(info_dict.get('ligand_mol'), size=(200, 200))
    _img.save(_imgdata, format='PNG')
    output_worksheet.merge_range('D2:E3', '', center_fmt)
    output_worksheet.insert_image('D2', 'f', {'image_data': _imgdata, 'x_offset': 35, 'y_offset': 7})

    output_worksheet.merge_range('F2:F3', 'Best', bold_fmt)
    output_worksheet.merge_range('G2:G3', info_dict.get('best_pdb'), center_fmt)
    
    cols = df.columns.tolist()
    output_worksheet.write(3, 0, 'Protein', subtitle_fmt)
    output_worksheet.write_row(3, 1, data=cols[1:], cell_format=subtitle_fmt)
    # cols.remove('PDB')
    # cols.remove('Ligand')
    # cols.remove('Structure')

    # output_worksheet.write(3, 0, 'Protein', subtitle_fmt)
    # output_worksheet.write(3, 1, 'PDB', subtitle_fmt)
    # output_worksheet.write(3, 2, 'Reference_Ligand', subtitle_fmt)
    # output_worksheet.write(3, 3, 'Structure', subtitle_fmt)
    # output_worksheet.write_row(3, 4, cols, subtitle_fmt)

    for i, (index, row) in enumerate(df.iterrows()):
        output_worksheet.set_row(i+4, 114, center_fmt)
        imgdata = BytesIO()

        output_worksheet.write(i+4, 0, row['Receptor'])
        output_worksheet.write(i+4, 1, row['PDB'])
        output_worksheet.write(i+4, 2, row['Reference_Ligand'])
        try:
            img = Draw.MolToImage(row['Structure'], size=(150, 150))
            img.save(imgdata, format='PNG')
            output_worksheet.insert_image(i+4, 3, "f", {'image_data': imgdata, 'y_offset': 1, 'x_offset': 1})
        except Exception as e:
            pass
        

        # for col in cols:
        #     try:
        #         output_worksheet.write(i+4, cols.index(col)+4, row[col], center_fmt)
        #     except Exception:
        #         pass
        try:
            output_worksheet.write(i+4, 4, row['Docking_Score'], center_fmt)
            output_worksheet.write(i+4, 5, row['rmsd'], center_fmt)
        except Exception:
            pass
    
        try:
            if row[cols[-1]] <= row['Docking_Score']:
                output_worksheet.write(i+4, 6, row[cols[-1]], highlight_fmt)
            else:
                output_worksheet.write(i+4, 6, row[cols[-1]], center_fmt)
        except Exception:
            pass
    
    for i in range(7):
        output_worksheet.write(len(df)+4, i, '', endline_fmt)
    output_workbook.close()

def _get_ligand_info_dict(ligand_file_path):

    if not os.path.exists(ligand_file_path):
        logger.warning(f'{ligand_file_path} not exists.')
    ligand_file_path = os.path.abspath(ligand_file_path)
    ligand_smiles = _get_smiles(ligand_file_path)
    ligand_mol = MolFromSmiles(ligand_smiles)

    info_dict = {}
    info_dict['ligand_name'] = ligand_file_path.split('/')[-1].split('.')[0]
    info_dict['ligand_smiles'] = ligand_smiles
    info_dict['ligand_mol'] = ligand_mol

    return info_dict

def generate_report(mappings:list, refefence_datalist:list, dock_datalist:list, ligand_save_dir:str, report_save_dir:str):
    '''
    为每个Ligand生成对接报告

    Parameters
    ----------
    mappings : list
        原始输入文件的映射关系
    refefence_datalist : list[dict]
        参考数据列表
    dock_datalist : list[dict]
        对接数据列表
    ligand_save_dir : str
        Ligand文件保存目录(用于读取结构)
    report_save_dir : str
        报告xlsx保存目录
    '''
    # 创建总报告模板
    template_df = pd.DataFrame(mappings)
    template_df.columns = ['Receptor', 'PDB', 'Ligand']
    template_df = template_df.drop_duplicates()
    template_df['Reference_Ligand'] = template_df['PDB'] + '-lig-' + template_df['Ligand']
    template_df = template_df[['Receptor', 'PDB', 'Reference_Ligand']]

    # 处理参考数据
    _redock_data = []
    _redock_df = pd.DataFrame(refefence_datalist)
    redock_df_group = _redock_df.groupby('PDB')
    for _pdb, pdb_df in redock_df_group:
        for index, row in pdb_df.iterrows():
            if row['Ligand'].startswith(row['PDB']):
                _redock_data.append(row[['PDB', 'Ligand', 'Docking_Score', 'rmsd']].to_dict())
    redock_df = pd.DataFrame(_redock_data)
    redock_df['Reference_Ligand'] = redock_df['Ligand']
    redock_df = redock_df[['PDB', 'Reference_Ligand', 'Docking_Score', 'rmsd']]
    
    final_df = pd.merge(template_df, redock_df, on=['PDB', 'Reference_Ligand'], how='left')
    final_df = _generate_img(final_df, 'Reference_Ligand', ligand_save_dir)

    # 处理对接数据
    dock_df = pd.DataFrame(dock_datalist)
    ligand_df_group = dock_df.groupby('Ligand')

    if len(ligand_df_group) == 0:
        logger.warning('No docking data found.')
        return
    elif len(ligand_df_group) > 10:
        logger.warning('Too many docking data found.')
        if not input('Continue?[y/n]').lower() == 'y':
            return

    for _ligand, ligand_df in ligand_df_group:
        
        info_dict = _get_ligand_info_dict(os.path.join(ligand_save_dir, f'{_ligand}.mae'))
        ligand_df['Reference_Ligand'] = ligand_df['PDB'] + '-lig-' + ligand_df['Original']
        ligand_df = ligand_df[['Reference_Ligand', 'Docking_Score']]
        ligand_df[_ligand] = ligand_df['Docking_Score']
        ligand_df = ligand_df.drop(columns=['Docking_Score'])
        _report_df = final_df.merge(ligand_df, how='left', on='Reference_Ligand')

        info_dict['higher_count'] = _report_df[_report_df[_ligand] <= _report_df['Docking_Score']].shape[0]
        info_dict['lower_count'] = _report_df[_report_df[_ligand] > _report_df['Docking_Score']].shape[0]
        info_dict['best_pdb'] = _report_df.loc[_report_df[_ligand].idxmin(), 'PDB']
        _write_xlsx(_report_df, os.path.join(report_save_dir, f'{_ligand}_report.xlsx'), info_dict)
