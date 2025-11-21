import json
import logging
import os
from urllib import parse

import pandas as pd
import requests
import yaml

from pyCADD.Demand.config import IGNORE_LIG, BaseQueryCfg, BaseQueryPDB, QueryConfigKeys
from pyCADD.utils.common import FixedConfig
from pyCADD.utils.tool import read_file, write_file

logger = logging.getLogger(__name__)

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/"
PDB_URL = "https://data.rcsb.org/graphql?query="
PDB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query?json="


def parse_uniport(uniprot_file_path: str):
    data = json.loads(read_file(uniprot_file_path))
    data = pd.DataFrame(data["uniProtKBCrossReferences"])
    try:
        pdb_list = data["id"].to_list()
    except KeyError:
        raise ValueError(
            f'Uniprot ID {os.path.basename(uniprot_file_path).split(".")[0]} has no PDB crystal data.'
        )
    return pdb_list


def query_uniprot_id_on_pdb(uniprot_id: str, save_dir: str = None) -> list:
    query = BaseQueryPDB(uniprot_id=uniprot_id).get_query()
    query = parse.quote(query)
    res = requests.get(PDB_SEARCH_URL + query)
    if res.status_code != 200:
        raise ValueError(f"Query failed: {res.text}")
    if save_dir is not None:
        write_file(os.path.join(save_dir, uniprot_id + ".json"), res.text)

    search_result = res.json()
    result_length = int(search_result["total_count"])
    pdb_list = [result["identifier"] for result in search_result["result_set"]]
    if len(pdb_list) != result_length:
        raise ValueError(
            f"Uniprot ID {uniprot_id} has {result_length} PDB crystal data, but only {len(pdb_list)} is returned."
        )
    return pdb_list


def query_pdb(pdb_list: list, save_path: str = None, query_cfg: BaseQueryCfg = None):
    query_cfg = BaseQueryCfg(pdb_list) if query_cfg is None else query_cfg(pdb_list)
    query = query_cfg.get_query()
    query = parse.quote(query)
    try:
        res = requests.get(PDB_URL + query)
    except requests.exceptions.SSLError:
        raise requests.exceptions.SSLError("SSL Error: Please check your network connection.")
    if res.status_code != 200:
        raise ValueError(f"Query failed: {res.text}")
    data_json = res.json()
    if save_path is not None:
        write_file(save_path, json.dumps(data_json))
    return data_json


def get_nested_value(data, keys):
    try:
        for key in keys:
            if isinstance(data, list):
                data = data[key]
            elif isinstance(data, dict):
                data = data.get(key, None)
            if data is None:
                return None
        return data
    except KeyError:
        return None


class QueryClient:
    def __init__(
        self, uniprot_id: str = None, pdb_list_file: str = None, save_dir: str = None
    ) -> None:
        self.uniprot_id = uniprot_id
        self.pdb_list_file = pdb_list_file
        self.save_dir = save_dir or os.getcwd()

        if self.uniprot_id is None and self.pdb_list_file is None:
            raise ValueError("Uniprot ID or PDB list file must be provided.")

        self.pdb_list = None
        self.pdb_data = None
        self.pdb_query_cfg = None
        self.data_dict = None
        self.mutation_pdb = None
        self.pairs_clean = None
        self.output_data = None
        self.apo = None
        self.from_file = False

        if self.pdb_list_file is not None:
            self._parse_pdb_file(self.pdb_list_file)
            self.uniprot_id = os.path.basename(self.pdb_list_file).split(".")[0]
            self.from_file = True

    @property
    def uniprot_save_dir(self):
        uniprot_save_dir = os.path.join(self.save_dir, "query_data", "uniprot")
        os.makedirs(uniprot_save_dir, exist_ok=True)
        return uniprot_save_dir

    @property
    def pdb_save_dir(self):
        pdb_save_dir = os.path.join(self.save_dir, "query_data", "pdb")
        os.makedirs(pdb_save_dir, exist_ok=True)
        return pdb_save_dir

    def _parse_pdb_file(self, pdb_file_path: str, pdbid_column: str = None):
        if not os.path.exists(pdb_file_path):
            raise FileNotFoundError(f"PDB file {pdb_file_path} not found.")
        pdbid_column = "PDBID" if pdbid_column is None else pdbid_column
        df = pd.read_csv(pdb_file_path)
        try:
            self.pdb_list = df[pdbid_column].to_list()
        except KeyError:
            raise KeyError(f"PDB file {pdb_file_path} has no column {pdbid_column}.")

    def get_apo(self):
        return self.apo

    def get_mutations(self):
        return self.mutation_pdb

    def _query_uniprot(self):
        save_path = os.path.join(self.uniprot_save_dir, self.uniprot_id + ".json")
        res = requests.get(UNIPROT_URL + self.uniprot_id + "?format=json&fields=xref_pdb")
        if res.status_code != 200:
            raise ValueError(f"Uniprot ID: {self.uniprot_id} is not found.")
        write_file(save_path, res.text)

    def _query_uniprot_id(self):
        self.pdb_list = query_uniprot_id_on_pdb(self.uniprot_id, save_dir=self.uniprot_save_dir)

    def _query_pdb(self):
        if self.pdb_list is None:
            self._query_uniprot_id()
        save_path = os.path.join(self.pdb_save_dir, self.uniprot_id + ".json")
        self.data_dict = query_pdb(self.pdb_list, save_path, self.pdb_query_cfg)
        if self.data_dict is None:
            raise ValueError(f"Uniprot ID: {self.uniprot_id} has no PDB crystal data.")
        self._parse_json()

    def _parse_json(self):
        parse_data = []
        js_data = get_nested_value(self.data_dict, QueryConfigKeys.DATA_KEYS)

        for dic_pdb in js_data:
            d = {}
            d["PDBID"] = get_nested_value(dic_pdb, QueryConfigKeys.PDBID_KEYS)  # PDB ID
            d["title"] = get_nested_value(dic_pdb, QueryConfigKeys.TITLE_KEYS)  # 标题
            d["resolution"] = get_nested_value(dic_pdb, QueryConfigKeys.RESOLUTION_KEYS)  # 分辨率
            d["reference"] = get_nested_value(dic_pdb, QueryConfigKeys.REFERENCE_KEYS)  # 参考文献
            d["authors"] = "".join(
                i + " " for i in get_nested_value(dic_pdb, QueryConfigKeys.AUTHOR_KEYS)
            )
            d["DOI"] = get_nested_value(dic_pdb, QueryConfigKeys.REF_DOI_KEYS)  # DOI
            d["polymer_entity"] = get_nested_value(
                dic_pdb, QueryConfigKeys.POLYMER_ENTITY_NUM_KEYS
            )  # 多肽体实体数
            d["nonpolymer_entity"] = get_nested_value(
                dic_pdb, QueryConfigKeys.NONPOLYMER_ENTITY_NUM_KEYS
            )  # 非多肽体实体数

            # 解析多肽体实体
            polymer_entities = get_nested_value(dic_pdb, QueryConfigKeys.POLYMER_ENTITY_KEYS)
            for index, polymer in enumerate(polymer_entities):
                d[f"polymer_entities_{index}_name"] = get_nested_value(
                    polymer, QueryConfigKeys.POLYMER_NAME_KEYS
                )  # 名称
                d[f"polymer_entities_{index}_type"] = get_nested_value(
                    polymer, QueryConfigKeys.POLYMER_TYPE_KEYS
                )  # 类型
                d[f"polymer_entities_{index}_mutation_num"] = get_nested_value(
                    polymer, QueryConfigKeys.POLYMER_MUTATION_NUM_KEYS
                )  # 突变数
                d[f"polymer_entities_{index}_mutation"] = get_nested_value(
                    polymer, QueryConfigKeys.POLYMER_MUTATION_KEYS
                )  # 突变数
                d[f"polymer_entities_{index}_chain"] = get_nested_value(
                    polymer, QueryConfigKeys.POLYMER_CHAIN_ID_KEYS
                )  # 所在链
                try:
                    d["polymer_entities_" + str(index) + "_source_organism"] = (
                        ",".join(
                            list(d.values())[0]
                            for d in polymer["rcsb_entity_source_organism"][0]["rcsb_gene_name"]
                        )
                        if polymer["rcsb_entity_source_organism"]
                        else None
                    )  # 结构来源蛋白
                except TypeError:
                    d["polymer_entities_" + str(index) + "_source_organism"] = None

            # 解析非多肽体实体
            nonpolymer_entities = get_nested_value(dic_pdb, QueryConfigKeys.NONPOLYMER_ENTITY_KEYS)
            if nonpolymer_entities is not None:
                for index, nonpolymer in enumerate(nonpolymer_entities):
                    d[f"nonpolymer_entities_{index}_name"] = get_nested_value(
                        nonpolymer, QueryConfigKeys.NONPOLYMER_NAME_KEYS
                    )  # 名称
                    d[f"nonpolymer_entities_{index}_ID"] = get_nested_value(
                        nonpolymer, QueryConfigKeys.NONPOLYMER_ID_KEYS
                    )  # ID
                    d[f"nonpolymer_entities_{index}_SMILES"] = get_nested_value(
                        nonpolymer, QueryConfigKeys.NONPOLYMER_SMILES_KEYS
                    )  # SMILES
                    d[f"nonpolymer_entities_{index}_chain"] = ",".join(
                        i
                        for i in get_nested_value(
                            nonpolymer, QueryConfigKeys.NONPOLYMER_CHAIN_ID_KEYS
                        )
                    )

            parse_data.append(d)
        self.pdb_data = parse_data

    def query(self):
        self._query_pdb()
        csv_file = f"{os.path.abspath(os.path.join(self.pdb_save_dir, self.uniprot_id))}.csv"
        pd.DataFrame(self.pdb_data).to_csv(csv_file, index=False)
        logger.info(f"Parsed PDB data is saved to {csv_file}")

    def clean_pdb_data(
        self, del_mutations: bool = True, del_ignore_lig: bool = True, cutoff: float = None
    ):
        """
        清洗 pdb 数据:
            * 去除Apo晶体
            * 去除配体未结合于目标链的晶体
            * 去除非WideType晶体(optional)
            * 去除非配体的小分子(e.g. DMS, optional)
            * 去除分辨率高于Cutoff的晶体(optional)
        """
        self.apo = []
        total_result = {}
        if self.from_file:
            raise ValueError("PDB list is from file, cannot be cleaned.")
        if self.data_dict is None:
            self.query()
        data_list = get_nested_value(self.data_dict, QueryConfigKeys.DATA_KEYS)

        for pdb in data_list:
            target_ligs = []
            pdbid = get_nested_value(pdb, QueryConfigKeys.PDBID_KEYS)
            if get_nested_value(pdb, QueryConfigKeys.NONPOLYMER_ENTITY_KEYS) is None:
                self.apo.append(pdbid)
                continue
            if cutoff is not None:
                resolution = get_nested_value(pdb, QueryConfigKeys.RESOLUTION_KEYS)
                if float(resolution) > cutoff:
                    continue
                elif resolution is None:
                    logger.warning(f"PDB {pdbid} has no resolution data.")
                    continue

            # 目标蛋白所在链
            for poly_entity in get_nested_value(pdb, QueryConfigKeys.POLYMER_ENTITY_KEYS):
                _uniprot_id = get_nested_value(poly_entity, QueryConfigKeys.POLYMER_UNIPROT_ID_KEYS)
                if isinstance(_uniprot_id, list):
                    if _uniprot_id[0] == self.uniprot_id:
                        target_chains = set(
                            get_nested_value(
                                poly_entity, QueryConfigKeys.POLYMER_CHAIN_ID_KEYS
                            ).split(",")
                        )

            for nonpoly_entity in get_nested_value(pdb, QueryConfigKeys.NONPOLYMER_ENTITY_KEYS):
                binding_chains = set(
                    get_nested_value(nonpoly_entity, QueryConfigKeys.NONPOLYMER_CHAIN_ID_KEYS)
                )
                # 存在并集 即结合在目标链上
                if binding_chains.intersection(target_chains):
                    target_ligs.append(
                        get_nested_value(nonpoly_entity, QueryConfigKeys.NONPOLYMER_ID_KEYS)
                    )
            total_result[pdbid] = target_ligs

        # Filter Apo
        self.pairs_clean = {k: v for k, v in total_result.items() if set(v).difference(IGNORE_LIG)}
        self.output_data = self.pairs_clean

        if del_mutations:
            mutations = self.get_mutation_pdb()
            self.output_data = {k: v for k, v in self.output_data.items() if k not in mutations}

        if del_ignore_lig:
            self.output_data = {
                key: [v for v in value if v not in IGNORE_LIG]
                for key, value in self.output_data.items()
            }
        return self.output_data

    def get_mutation_pdb(self):
        mutations = []
        data_dict = get_nested_value(self.data_dict, QueryConfigKeys.DATA_KEYS)
        for pdb in data_dict:
            pdbid = get_nested_value(pdb, QueryConfigKeys.PDBID_KEYS)
            for entity in get_nested_value(pdb, QueryConfigKeys.POLYMER_ENTITY_KEYS):
                if get_nested_value(entity, QueryConfigKeys.POLYMER_TYPE_KEYS) != "Protein":
                    continue
                if get_nested_value(entity, QueryConfigKeys.POLYMER_MUTATION_NUM_KEYS) != 0:
                    mutations.append(pdbid)
        self.mutation_pdb = set(mutations)
        return self.mutation_pdb

    def generate_inputfile(self, path: str, _format: str = None):
        if _format is None:
            ext = path.split(".")[-1].lower()
            if ext in ["csv", "txt"]:
                _format = "csv"
            elif ext in ["in", "ini"]:
                _format = "ini"
            elif ext in ["yml", "yaml"]:
                _format = "yaml"
            else:
                raise NotImplementedError(f'Unsupported format: {path.split(".")[-1]}')

        output_data = self.output_data.copy()
        if _format == "yaml":
            for k, v in output_data.items():
                if len(v) == 1:
                    output_data[k] = v[0]
            with open(path, "w") as f:
                yaml.dump({self.uniprot_id: output_data}, f)

        elif _format == "csv":
            with open(path, "w") as f:
                for pdb, ligs in output_data.items():
                    for lig in ligs:
                        f.write(f"{pdb},{lig}\n")

        elif _format == "ini" or _format == "in":
            config = FixedConfig()
            config.add_section(self.uniprot_id)
            for pdb, ligs in output_data.items():
                ligs = ",".join(ligs)
                config.set(self.uniprot_id, pdb, ligs)
            with open(path, "w") as f:
                config.write(f)

        logger.info(f"Dock input file is saved to {path}.")
