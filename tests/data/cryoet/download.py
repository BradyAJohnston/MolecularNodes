import argparse

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

parser = argparse.ArgumentParser()
parser.add_argument(
    "--download-segmask", action="store_true", help="Download segmentation mask files"
)
args = parser.parse_args()

client = Client()
run = Run.get_by_id(client, 506)

for annotation in run.annotations:
    # Always download OrientedPoint
    annotation.download(shape="OrientedPoint")
    # Optionally download SegmentationMask
    if args.download_segmask:
        annotation.download(shape="SegmentationMask")
