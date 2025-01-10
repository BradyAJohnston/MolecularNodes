# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "cryoet-data-portal",
# ]
# ///
from cryoet_data_portal import (
    Client,
    Run,
)

client = Client()
run = Run.get_by_id(client, 506)

for annotation in run.annotations:
    # Download only ndjson and mrc files
    for shape in ["OrientedPoint"]:
        annotation.download(shape=shape)
