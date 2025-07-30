import multiprocessing
import os
import re
from pathlib import Path

import requests

_FIES_URL = "https://kirnis.kapsi.fi/NOT2023/data/"
_FITS_RE = re.compile(r'<a href="(.*.fits)">.*</a>')
_FOLDERS_RE = re.compile(r'<a href="(.*/)">.*/</a>')


def _crawl_fits_links(url: str) -> list[str]:
    page_text = requests.get(url).text
    links = []

    for link in _FITS_RE.findall(page_text):
        links.append(url + link)

    for folder in _FOLDERS_RE.findall(page_text):
        sub_links = _crawl_fits_links(url + folder)
        links.extend(sub_links)

    return links


def _download_file(base_url: str, base_path: Path, file_url: str):
    file_path = base_path / (file_url.removeprefix(base_url))
    file_path.parent.mkdir(parents=True, exist_ok=True)

    response = requests.get(file_url, stream=True)
    if response.status_code != 200:
        return None

    with open(file_path, mode="wb") as file:
        for chunk in response.iter_content(chunk_size=2**16):
            file.write(chunk)
    return file_url


def _downloader_task(data):
    base_url, base_path, file_url = data
    return _download_file(base_path=base_path, base_url=base_url, file_url=file_url)


def _download_files(links: list[str], base_url: str, base_path: Path) -> tuple[set, set]:
    arguments = [(base_url, base_path, link) for link in links]
    succeeded = set()
    with multiprocessing.Pool(processes=os.cpu_count()) as pool:
        for url in pool.imap(_downloader_task, arguments):
            if url:
                print("\t\t * Downloaded:", url)
                succeeded.add(url)
    failed = set(links).difference(succeeded)
    return succeeded, failed


def download_spectra(out_folder: Path):
    print("### Downloading FIES data ###")

    print("\t * Crawling url", _FIES_URL)
    fies_urls = _crawl_fits_links(_FIES_URL)
    print(f"\t * Found {len(fies_urls)} urls")

    print(f"\t * Downloading files into {out_folder}")
    succeeded, failed = _download_files(links=fies_urls, base_url=_FIES_URL, base_path=out_folder)
    print(f"\t * {len(succeeded)} succeeded, {len(failed)} failed")
    for f in failed:
        print(f"\t\t - Failed: {f}")


if __name__ == "__main__":
    _out_folder = Path.cwd() / "data/fits/"
    download_spectra(_out_folder)
