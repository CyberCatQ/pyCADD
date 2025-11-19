FORCE_FILED_DICT = {
    "OPLS3e": "3",
    "OPLS4": "S-OPLS",
    "OPLS2005": "2005",
    "OPLS_2005": "2005",
    "OPLS3": "2.1",
}


class DefaultDataConfig:
    def __init__(self) -> None:
        self.properties: dict = {
            "Docking_Score": "r_i_docking_score",
            "rmsd": "r_i_glide_rmsd_to_input",
            "ligand_efficiency": "r_i_glide_ligand_efficiency",
            "rotatable_bonds": "i_i_glide_rotatable_bonds",
            "ecoul": "r_i_glide_ecoul",
            "evdw": "r_i_glide_evdw",
            "emodel": "r_i_glide_emodel",
            "energy": "r_i_glide_energy",
            "einternal": "r_i_glide_einternal",
        }

    def items(self):
        return self.properties.items()

    def keys(self):
        return self.properties.keys()

    def values(self):
        return self.properties.values()

    def get(self, key):
        return self.properties.get(key)

    def __getitem__(self, key):
        return self.properties[key]

    def __setitem__(self, key, value):
        return self.properties.__setitem__(key, value)


class SPConfig(DefaultDataConfig):
    def __init__(self) -> None:
        super().__init__()
        self.properties.update(
            {
                "lipo": "r_i_glide_lipo",
                "hbond": "r_i_glide_hbond",
                "metal": "r_i_glide_metal",
                "rewards": "r_i_glide_rewards",
                "erotb": "r_i_glide_erotb",
                "esite": "r_i_glide_esite",
            }
        )


class XPConfig(DefaultDataConfig):
    def __init__(self) -> None:
        super().__init__()
        self.properties.update({"XP_Hbond": "r_glide_XP_HBond"})


class DataConfig(DefaultDataConfig):
    def __init__(self, precision: str = None, properties: dict = None) -> None:
        if not precision:
            config_class = DefaultDataConfig
        elif precision.upper() == "SP":
            config_class = SPConfig
        elif precision.upper() == "XP":
            config_class = XPConfig

        self.properties = config_class().properties
        if isinstance(properties, dict):
            self.properties.update(properties)
