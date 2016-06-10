import os
import sys
import platform
import shutil

import requests


REQUIRED_SOFTWARE = ('bigWigAverageOverBed', 'validateFiles', )


def get_bin_root():
    root_dl = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'bin')
    if not os.path.exists(root_dl):
        os.makedirs(root_dl)
    return root_dl


def binaries_exist():
    root_dl = get_bin_root()
    return all([
        os.path.exists(os.path.join(root_dl, fn))
        for fn in REQUIRED_SOFTWARE
    ])


def get_root_url():
    system = platform.system()
    if system == 'Darwin':
        return 'http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64'
    elif system == 'Linux':
        return 'http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64'
    else:
        raise OSError('Mac or Linux system required.')


def download_ucsc_tools():
    bits, _ = platform.architecture()

    # check architecture
    if bits != '64bit':
        raise OSError('64-bit architecture required.')

    root_url = get_root_url()
    root_dl = get_bin_root()

    # download required software and place in appropriate location
    sys.stdout.write("Downloading UCSC binaries for ORIO\n")
    for fn in REQUIRED_SOFTWARE:
        url = os.path.join(root_url, fn)
        path = os.path.join(root_dl, fn)

        sys.stdout.write("Downloading: {}\n".format(url))
        r = requests.get(url, stream=True)
        if r.status_code == 200:
            with open(path, 'wb') as f:
                r.raw.decode_content = True
                shutil.copyfileobj(r.raw, f)
                os.chmod(path, 0o751)
        else:
            sys.stderr.write("URL returned a non-200 status\n")

    sys.stdout.write("Downloads complete!\n")

if __name__ == '__main__':
    download_ucsc_tools()
