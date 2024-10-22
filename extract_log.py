import re
import glob
from typing import List
from pathlib import Path
import json


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


class PlatformLog:
    def __init__(self, path: Path):
        self.path = path
        self.urls = extract_download_urls(self.path)
        self.url_dict = self.urls_to_dict(self.urls)

    def urls_to_dict(self, urls: List[str]) -> dict:
        return {url.split("/")[-1]: url for url in urls}

    @property
    def platform(self) -> str:
        return self.path.name.removesuffix("_pip_log.txt")


class LogHandler:
    def __init__(self, logs: List[str]):
        self.log_filenames = [Path(log) for log in logs]
        self.log_dict: dict = {}

    def process_logs(self):
        for path in self.log_filenames:
            log = PlatformLog(path)
            self.log_dict[log.platform] = log.url_dict


handler = LogHandler(logs)
handler.process_logs()

with open("download_urls.json", "w") as f:
    json.dump(handler.log_dict, f, indent=True)
