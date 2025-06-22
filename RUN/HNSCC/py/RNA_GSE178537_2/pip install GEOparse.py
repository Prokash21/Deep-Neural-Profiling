import GEOparse

import GEOparse
gse = GEOparse.get_GEO(geo="GSE178537", destdir=".")


for gsm_name, gsm in gse.gsms.items():
    print(f"{gsm_name}: {gsm.metadata['title'][0]}")

import pandas as pd

# Create a DataFrame of GSM metadata
sample_data = pd.DataFrame({gsm_name: gsm.metadata for gsm_name, gsm in gse.gsms.items()}).T

# Display sample_data with selected columns (like title or characteristics)
print(sample_data[['title', 'source_name_ch1']])
