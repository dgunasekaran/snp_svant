import glob
import pandas as pd
from snakemake.io import expand
from snakemake.remote import FTP
from snakemake.utils import validate

ftp = FTP.RemoteProvider()

validate(config, schema="../schemas/config.schemas.yaml")

input_fp = config["samples"]["metadata"]
pd_header = 0 if config["samples"]["header"] else None


def get_samples():
    input_pd = pd.read_csv(input_fp, sep="\t", header=pd_header, names=["sample_name"])
    return input_pd

SAMPLES = get_samples()["sample_name"]
wildcard_constraints:
    sample="|".join(get_samples()["sample_name"])
