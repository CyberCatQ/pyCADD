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
			pdbx_vrpt_summary{
				PDB_resolution
			}
			rcsb_entry_info {
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