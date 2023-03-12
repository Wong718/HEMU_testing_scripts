"""
TFHeatmapSampleConfiguration.py
Stores sample information used for the creation of TF heatmap
"""
zea_dict = {
    "root": ["SRR13178424", "SRR13518146", "SRR13518151"],
    "ear": ["SRR12616469", "SRR7970747", "SRR7970750"],
    "leaf": ["SRR13791451", "SRR13791456", "SRR13791460"],
    "shoot": ["SRR2078285", "SRR2078286", "SRR2078289"],
    "kernel": ["SRR13077087", "SRR13077088", "SRR13077089", "SRR13077090"],
    "grain": ["SRR13002214", "SRR13002217", "SRR13002219"],
    "anther": ["SRR10855045", "SRR10855047", "SRR10855049"],
    "seed": ["ERR4691103", "ERR4691105", "ERR4691109"],
    "pericarp": ["SRR13664378", "SRR13664381", "SRR13664384"],
    "sheath": ["SRR7469475", "SRR7469478", "SRR7469488"],
    "blade": ["SRR7469496", "SRR7469498", "SRR7469500"],
    "ligule": ["SRR7469501", "SRR7469502", "SRR7469506"],
    "embryo": ["SRR513472", "SRR513468", "SRR513466"],
    "endosperm": ["SRR513461", "SRR513463", "SRR513464"],
}

coix_dict = {
    "leaf": ["SRR10193329", "SRR10193335", "SRR10193338"],
    "seed": ["SRR10208252", "SRR10208254", "SRR10208265"],
    "root": ["SRR10208255", "SRR10208256", "SRR10208257"],
    "glume": ["SRR9112527", "SRR9112528", "SRR9112529"],
    "male_flower": ["SRR9112530", "SRR9112531", "SRR9112543"],
}


def return_sample(species):
    if species == "coix":
        return coix_dict
    elif species == "zea":
        return zea_dict
    else:
        return None