# 非配体的小分子
IGNORE_LIG = ['EDO', 'DMS', 'IPA', 'TBY', 'ARS', 'EU', 'MG', 'IOD', 'ACT', 'CA', 'CAC', 'K', 'FMT', 'BU3', 'PGO', 'PE4',
               'PO4', 'BR', 'NO3', 'BCT', 'ZN', 'SO4', 'CL', 'NA', 'AU', 'GOL', 'NI', 'YT3', 'PEG', 'PGE']

class BaseQueryCfg:
	def __init__(self, pdb_list):
		pdb_list = ['"' + pdb + '"' for pdb in pdb_list]
		self.pdbs = ','.join(pdb_list)
		self.query = self._default_query
    
	def get_query(self):
		return self.query
	
	@property
	def _default_query(self):
		return '''{
		entries(entry_ids: [%s])
		{
			rcsb_id
			rcsb_primary_citation {
				pdbx_database_id_DOI
				rcsb_authors
				title
			}
			rcsb_binding_affinity {
				comp_id
				value
				provenance_code
				type
				unit
			}
			rcsb_entry_info {
       			resolution_combined
				deposited_nonpolymer_entity_instance_count
				deposited_polymer_entity_instance_count
			}
			struct {
				title
			}
			polymer_entities {
				entity_poly {
					rcsb_entity_polymer_type
					rcsb_mutation_count
					pdbx_strand_id
				}
				rcsb_polymer_entity_container_identifiers{
					uniprot_ids
				}
				rcsb_entity_source_organism {
					rcsb_gene_name {
						value
					}
				}
				rcsb_polymer_entity {
					pdbx_description
					pdbx_mutation
				}
			}
			nonpolymer_entities {
				nonpolymer_comp {
					chem_comp {
						id
						name
					}
					rcsb_chem_comp_descriptor {
						SMILES
					}
				}
				rcsb_nonpolymer_entity_container_identifiers {
					auth_asym_ids
						}
					}
				}
		}
		''' % self.pdbs

class BaseQueryPDB:
	def __init__(self, uniprot_id:str) -> None:
		self.uniprot_id = uniprot_id.strip()
  
	def get_query(self):
		return self._search_pdb_list

	@property
	def _search_pdb_list(self):
		return '''{
		"query": {
			"type": "group",
			"nodes": [
			{
				"type": "terminal",
				"service": "text",
				"parameters": {
				"attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
				"operator": "in",
				"value": [
					"%s"
				]
				}
			},
			{
				"type": "terminal",
				"service": "text",
				"parameters": {
				"attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
				"operator": "exact_match",
				"value": "UniProt"
				}
			}
			],
			"logical_operator": "and",
			"label": "nested-attribute"
		},
		"return_type": "entry",
		"request_options": {
		"paginate": {
			"start": 0,
			"rows": 10000
				}
			}
		}''' % self.uniprot_id
  
DATA_KEYS = ['data', 'entries']
PDBID_KEYS = ['rcsb_id']
TITLE_KEYS = ['struct', 'title']
RESOLUTION_KEYS = ['rcsb_entry_info', 'resolution_combined', 0]
AUTHOR_KEYS = ['rcsb_primary_citation', 'rcsb_authors']
REFERENCE_KEYS = ['rcsb_primary_citation', 'title']
REF_DOI_KEYS = ['rcsb_primary_citation', 'pdbx_database_id_DOI']
POLYMER_ENTITY_NUM_KEYS = ['rcsb_entry_info', 'deposited_polymer_entity_instance_count']
NONPOLYMER_ENTITY_NUM_KEYS = ['rcsb_entry_info', 'deposited_nonpolymer_entity_instance_count']

POLYMER_ENTITY_KEYS = ['polymer_entities']
POLYMER_NAME_KEYS = ['rcsb_polymer_entity', 'pdbx_description']
POLYMER_TYPE_KEYS = ['entity_poly', 'rcsb_entity_polymer_type']
POLYMER_MUTATION_NUM_KEYS = ['entity_poly', 'rcsb_mutation_count']
POLYMER_MUTATION_KEYS = ['rcsb_polymer_entity', 'pdbx_mutation']
POLYMER_CHAIN_ID_KEYS = ['entity_poly', 'pdbx_strand_id']
POLYMER_UNIPROT_ID_KEYS = ['rcsb_polymer_entity_container_identifiers', 'uniprot_ids']

NONPOLYMER_ENTITY_KEYS = ['nonpolymer_entities']
NONPOLYMER_NAME_KEYS = ['nonpolymer_comp', 'chem_comp', 'name']
NONPOLYMER_ID_KEYS = ['nonpolymer_comp', 'chem_comp', 'id']
NONPOLYMER_SMILES_KEYS = ['nonpolymer_comp', 'rcsb_chem_comp_descriptor', 'SMILES']
NONPOLYMER_CHAIN_ID_KEYS = ['rcsb_nonpolymer_entity_container_identifiers', 'auth_asym_ids']
