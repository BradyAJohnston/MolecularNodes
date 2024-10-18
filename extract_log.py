import re
import glob


def extract_download_urls(log_file_path):
    pattern = r"Downloading link ([^ ]+)"
    download_urls = []

    try:
        with open(log_file_path, "r") as file:
            for line in file:
                match = re.search(pattern, line)
                if match:
                    download_urls.append(match.group(1))
    except Exception as e:
        print(f"Error reading {log_file_path}: {e}")

    return download_urls


logs = glob.glob("*log.txt")

all_urls = []
for log in logs:
    all_urls.extend(extract_download_urls(log))
all_urls.sort()

filenames = [url.split("/")[-1] for url in all_urls]

print(*zip(filenames, all_urls))
