import re
import tifffile
import glob

dir = r'/Users/eliasouellet-oviedo/Documents/Stage 2/Codes/Données/2024-06-18_TestSoftware/test save online-001/Bliq VMS/C1'  # Replace with your directory
files = glob.glob(f'{dir}/*.tif', recursive=True)
file = files[0]

print(tifffile.TiffFile(file).ome_metadata)

with tifffile.TiffFile(file) as tif:
    ome_metadata = tif.ome_metadata


def extract_value(tag: str, string: str) -> str | None:
    pattern = re.compile(rf'{tag}(\d+\.?\d*)')
    match = pattern.search(string)
    if match:
        return match.group(1)
    return None


def extract_float_value(tag: str, string: str) -> float | None:
    value = extract_value(tag, string)
    if value:
        try:
            return float(value)
        except ValueError:
            print(f"Could not convert {value} from {tag} to float.")
            return
    return


frame_rate = extract_float_value('<M K="FrameRate">', ome_metadata)
projection_count = extract_float_value('<M K="Projection count">', ome_metadata)

print(f"{frame_rate=}")
print(f"{projection_count=}")
