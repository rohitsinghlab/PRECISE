#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shutil
import stat
import sys
import tarfile
import tempfile
import zipfile
from pathlib import Path
from typing import Optional
from urllib.request import urlopen, Request

DEFAULT_APBS_URL = "https://github.com/Electrostatics/apbs/releases/download/v3.4.1/APBS-3.4.1.Linux.zip"
DEFAULT_MSMS_URL = "https://ccsb.scripps.edu/msms/download/933/"


def get_project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def set_executable(file_path: Path) -> None:
    mode = file_path.stat().st_mode
    file_path.chmod(mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def download(url: str, dest: Path, force: bool = False) -> None:
    if dest.exists() and not force:
        return
    ensure_dir(dest.parent)
    req = Request(url, headers={"User-Agent": "curl/7.81.0"})
    with urlopen(req) as response, tempfile.NamedTemporaryFile(delete=False) as tmp:
        shutil.copyfileobj(response, tmp)
        tmp_path = Path(tmp.name)
    dest.write_bytes(tmp_path.read_bytes())
    tmp_path.unlink(missing_ok=True)


def extract_zip(zip_file: Path, dest_dir: Path, force: bool) -> None:
    if dest_dir.exists() and force:
        shutil.rmtree(dest_dir)
    ensure_dir(dest_dir)
    with zipfile.ZipFile(zip_file, "r") as zf:
        zf.extractall(dest_dir)


def extract_tar_gz(tar_file: Path, dest_dir: Path, force: bool) -> None:
    if dest_dir.exists() and force:
        shutil.rmtree(dest_dir)
    ensure_dir(dest_dir)
    with tarfile.open(tar_file, "r:gz") as tf:
        tf.extractall(dest_dir)


def link_or_copy(source: Path, target: Path, force: bool) -> None:
    if target.exists() or target.is_symlink():
        if not force:
            return
        if target.is_dir() and not target.is_symlink():
            shutil.rmtree(target)
        else:
            target.unlink(missing_ok=True)
    try:
        target.symlink_to(source)
    except OSError:
        shutil.copy2(source, target) if source.is_file() else shutil.copytree(
            source, target
        )


def find_binary(root: Path, name: str) -> Optional[Path]:
    candidates = [p for p in root.rglob(f"{name}*") if p.is_file()]
    for p in candidates:
        if p.name == name:
            return p
    return sorted(candidates)[0] if candidates else None


def install_msms(install_dir: Path, url: str, force: bool) -> list[Path]:
    print(f"[MSMS] Using URL: {url}")
    packages_dir = install_dir / "packages" / "msms"
    bin_dir = install_dir / "bin"
    ensure_dir(packages_dir)
    ensure_dir(bin_dir)

    msms_binary = bin_dir / "msms"
    if msms_binary.exists() and os.access(str(msms_binary), os.X_OK) and not force:
        print("[MSMS] Already installed. Use --force to reinstall.")
        return [msms_binary]

    # The URL points to a tar.gz file
    archive = packages_dir / "msms.tar.gz"
    print(f"[MSMS] Downloading to {archive}...")
    download(url, archive, force)

    extract_root = packages_dir / "extracted"
    print(f"[MSMS] Extracting to {extract_root}...")
    extract_tar_gz(archive, extract_root, force)

    linked = []

    msms_path = None
    for item in extract_root.rglob("*"):
        if item.is_file() and "msms" in item.name.lower() and "x86_64" in item.name:
            msms_path = item
            break

    if not msms_path:
        msms_path = find_binary(extract_root, "msms")

    if not msms_path:
        raise RuntimeError(
            f"[MSMS] Could not locate MSMS binary in extracted content at {extract_root}"
        )

    print(f"[MSMS] Found MSMS binary: {msms_path}")

    # Link/copy as 'msms' (without the version suffix)
    target = bin_dir / "msms"
    link_or_copy(msms_path, target, force)
    set_executable(target)
    linked.append(target)

    print(f"[MSMS] Linked binary: msms")
    return linked


def install_apbs(install_dir: Path, url: str, force: bool) -> list[Path]:
    print(f"[APBS] Using URL: {url}")
    packages_dir = install_dir / "packages" / "apbs"
    bin_dir = install_dir / "bin"
    ensure_dir(packages_dir)
    ensure_dir(bin_dir)

    apbs_binary = bin_dir / "apbs"
    if apbs_binary.exists() and os.access(str(apbs_binary), os.X_OK) and not force:
        print("[APBS] Already installed. Use --force to reinstall.")
        return [apbs_binary]

    archive = packages_dir / (url.split("/")[-1] if "/" in url else "apbs.zip")
    print(f"[APBS] Downloading to {archive}...")
    download(url, archive, force)

    extract_root = packages_dir / "extracted"
    print(f"[APBS] Extracting to {extract_root}...")
    extract_zip(archive, extract_root, force)

    linked = []

    # Find bin directory and link its contents
    bin_dir_src = next((p for p in extract_root.rglob("bin") if p.is_dir()), None)
    if bin_dir_src:
        print(f"[APBS] Found bin directory at: {bin_dir_src}")
        for item in bin_dir_src.iterdir():
            if item.is_file():
                target = bin_dir / item.name
                link_or_copy(item, target, force)
                set_executable(target)
                linked.append(target)
    else:
        # No bin dir found, search for apbs directly
        apbs_path = find_binary(extract_root, "apbs")
        if not apbs_path:
            raise RuntimeError(
                "[APBS] Could not locate 'bin' or 'apbs' in extracted content."
            )
        target = bin_dir / "apbs"
        link_or_copy(apbs_path, target, force)
        set_executable(target)
        linked.append(target)

    # Search for multivalue if not already linked
    if not any(p.name == "multivalue" for p in linked):
        print(f"Searching for multivalue in {extract_root}")
        path = find_binary(extract_root, "multivalue")
        print(f"Found {path} for multivalue")
        if path:
            target = bin_dir / "multivalue"
            link_or_copy(path, target, force)
            set_executable(target)
            linked.append(target)

    print(f"[APBS] Linked binaries: {[p.name for p in linked]}")
    return linked


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Install APBS (with multivalue) and MSMS"
    )

    # Determine default installation directory
    # We default to the user directory to ensure tools are found regardless of
    # where the package is installed (e.g. source vs site-packages)
    default_dir = Path.home() / ".precise" / "tools"

    parser.add_argument(
        "--install-dir",
        type=Path,
        default=default_dir,
        help=f"Installation directory (default: {default_dir})",
    )
    parser.add_argument("--apbs-url", default=DEFAULT_APBS_URL)
    parser.add_argument("--msms-url", default=DEFAULT_MSMS_URL)
    parser.add_argument("--force", action="store_true", help="Force reinstall")
    parser.add_argument(
        "--skip-msms", action="store_true", help="Skip MSMS installation"
    )
    parser.add_argument(
        "--skip-apbs", action="store_true", help="Skip APBS installation"
    )
    args = parser.parse_args(argv)

    install_dir = args.install_dir.resolve()

    # Create directories if they don't exist
    try:
        ensure_dir(install_dir / "bin")
        ensure_dir(install_dir / "packages")
    except PermissionError:
        print(f"[ERROR] Permission denied: Cannot create directories in {install_dir}")
        return 1

    print(f"[INFO] Install directory: {install_dir}")

    all_linked = []

    if not args.skip_apbs:
        apbs_linked = install_apbs(install_dir, args.apbs_url, args.force)
        all_linked.extend(apbs_linked)

    if not args.skip_msms:
        msms_linked = install_msms(install_dir, args.msms_url, args.force)
        all_linked.extend(msms_linked)

    bin_dir = install_dir / "bin"
    print(f"\n[INFO] Installation complete. Binaries in: {bin_dir}")
    print(f'[INFO] Add to PATH: export PATH="{bin_dir}:$PATH"')

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
