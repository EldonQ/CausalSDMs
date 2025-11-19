#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
GDW（Global Dam Watch）数据自动抓取与下载脚本

功能概述（中文注释，便于本地二次开发）：
1. 从 GDW 数据库主页抓取页面内容与所有可见超链接；
2. 识别潜在可下载资源链接（常见数据后缀与 Google Drive 链接）；
3. 对网页型链接可按设定的最大深度继续抓取一层或多层，发现更多下载端点；
4. 对直接文件链接使用 requests 流式断点式下载，对 Google Drive 链接使用 gdown 下载；
5. 跳过已存在且非空的文件，支持失败重试与超时设置；
6. 输出完整的下载清单（manifest_gdw.csv）与日志，便于复现实验流程；

使用方法（Windows）：
    配置好 Python 环境与依赖后，直接在命令行执行：
        python gdw_download.py

注意：
    - 本脚本仅抓取公开可见链接；如需访问受限数据，请遵循原网站的许可与条款；
    - 若页面通过动态脚本加载或需要交互登录，可能需要手动补充下载端点或使用浏览器导出的链接；
    - 建议在稳定网络环境下运行，避免中途断线导致的下载失败。
"""

from __future__ import annotations

import csv
import datetime as dt
import logging
import os
import re
import sys
import time
from pathlib import Path
from typing import Iterable, List, Tuple, Dict, Set

# -------- 依赖检查 --------
try:
    import requests
    from bs4 import BeautifulSoup
    from tqdm import tqdm
    # lxml 是 beautifulsoup 高效解析所必需的
    import lxml
except ImportError as e:
    missing_module = str(e).split("'")[-2]
    print(f"错误: 必需的 Python 库 '{missing_module}' 未安装。", file=sys.stderr)
    print("请在命令行中运行以下命令来安装依赖:", file=sys.stderr)
    print("pip install requests beautifulsoup4 lxml tqdm gdown", file=sys.stderr)
    sys.exit(1)


# gdown 用于 Google Drive 链接的稳定下载
try:
    import gdown  # type: ignore
except Exception:  # noqa: E722
    gdown = None


# ----------------------------- 常量与全局配置 -----------------------------

# GDW 数据库主页（如站点结构更新，请据实际情况调整）
GDW_INDEX_URL = "https://www.globaldamwatch.org/database"

# 常见可下载数据的后缀集合（可按需扩展）
FILE_EXTENSIONS = (
    ".zip",
    ".7z",
    ".rar",
    ".csv",
    ".xlsx",
    ".xls",
    ".geojson",
    ".json",
    ".gpkg",
    ".tif",
    ".tiff",
    ".kmz",
    ".kml",
    ".parquet",
)

# 请求头，模拟常见浏览器，减少被误判为爬虫的概率
DEFAULT_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/124.0.0.0 Safari/537.36"
    ),
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    "Accept-Language": "en-US,en;q=0.9",
}


# ----------------------------- 工具函数：路径与日志 -----------------------------

def ensure_dir(path: Path) -> None:
    """确保目录存在，不存在则创建。"""
    path.mkdir(parents=True, exist_ok=True)


def setup_logging(log_dir: Path) -> None:
    """配置日志输出（文件 + 控制台）。"""
    ensure_dir(log_dir)
    log_file = log_dir / f"gdw_download_{dt.datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[
            logging.FileHandler(log_file, encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    logging.info("日志初始化完成：%s", log_file)


# ----------------------------- 网络与解析 -----------------------------

def http_get(url: str, timeout: int = 30) -> requests.Response:
    """发起 GET 请求，返回 Response；异常由上层处理。

    为了稳健性：
    - 使用默认请求头；
    - 由调用方控制重试；
    - 不在此处直接抛弃非 200 状态，交由上层判定。
    """
    session = requests.Session()
    session.headers.update(DEFAULT_HEADERS)
    resp = session.get(url, timeout=timeout, allow_redirects=True)
    return resp


def is_google_drive_url(url: str) -> bool:
    """判断是否为 Google Drive 链接（gdown 可处置）。"""
    return any(host in url for host in [
        "drive.google.com",
        "docs.google.com/uc",
        "googleusercontent.com",
    ])


def is_file_like_url(url: str) -> bool:
    """判断 URL 是否像是直接文件下载端点（根据后缀粗略判断）。"""
    lower = url.lower().split("?")[0].split("#")[0]
    return lower.endswith(FILE_EXTENSIONS)


def normalize_and_filter_links(base_url: str, links: Iterable[str]) -> List[str]:
    """标准化过滤链接：去空、去重复、去非 http(s)。"""
    from urllib.parse import urljoin

    normalized: List[str] = []
    seen: Set[str] = set()
    for href in links:
        if not href:
            continue
        href = href.strip()
        if href.startswith("#"):
            continue
        if href.startswith("/"):
            href = urljoin(base_url, href)
        if href.startswith("http://") or href.startswith("https://"):
            if href not in seen:
                seen.add(href)
                normalized.append(href)
    return normalized


def extract_links_from_html(html: str, base_url: str) -> List[str]:
    """从 HTML 中抽取所有 a 标签链接。"""
    soup = BeautifulSoup(html, "lxml")
    anchors = soup.find_all("a")
    hrefs: List[str] = []
    for a in anchors:
        try:
            href = a.get("href")
        except Exception:  # noqa: E722
            href = None
        if href:
            hrefs.append(href)
    return normalize_and_filter_links(base_url, hrefs)


def crawl_and_collect(
    seed_urls: List[str],
    max_depth: int,
    timeout: int,
) -> Tuple[List[str], List[str], List[str]]:
    """广度优先抓取链接，收集三类 URL：
    - direct_files: 直接可下载的文件端点（后缀匹配）
    - gdrive_files: Google Drive 链接（用 gdown 下载）
    - visited_pages: 实际访问过的网页（便于记录）

    max_depth 表示继续深入的网页层数（0 表示只用种子页）。
    """
    from collections import deque

    queue = deque([(u, 0) for u in seed_urls])
    visited_pages: Set[str] = set()
    direct_files: Set[str] = set()
    gdrive_files: Set[str] = set()

    while queue:
        url, depth = queue.popleft()
        if url in visited_pages:
            continue
        
        logging.info(f"正在分析页面 [深度 {depth}]: {url}")
        visited_pages.add(url)

        try:
            resp = http_get(url, timeout=timeout)
        except Exception as e:  # noqa: E722
            logging.warning("访问失败，将跳过：%s | 错误：%s", url, e)
            continue

        ctype = resp.headers.get("Content-Type", "").lower()
        if resp.status_code != 200:
            logging.warning("HTTP状态异常，将跳过：%s | 状态码：%s", url, resp.status_code)
            continue

        # 若返回的是二进制文件而非 HTML，直接判定为下载端点
        if ("text/html" not in ctype) and ("xml" not in ctype):
            if is_google_drive_url(url):
                gdrive_files.add(url)
            else:
                direct_files.add(url)
            continue

        # 解析 HTML 抽取子链接
        try:
            html = resp.text
        except Exception:  # noqa: E722
            try:
                html = resp.content.decode("utf-8", errors="ignore")
            except Exception:  # noqa: E722
                logging.warning("无法解析 HTML，将跳过：%s", url)
                continue

        child_links = extract_links_from_html(html, url)
        for lk in child_links:
            if is_google_drive_url(lk):
                gdrive_files.add(lk)
            elif is_file_like_url(lk):
                direct_files.add(lk)
            else:
                # 在深度限制内继续抓取子页面
                if depth < max_depth:
                    queue.append((lk, depth + 1))

    return sorted(direct_files), sorted(gdrive_files), sorted(visited_pages)


# ----------------------------- 下载实现 -----------------------------

def infer_filename_from_url(url: str) -> str:
    """根据 URL 推断文件名（保留查询字符串前的最后路径段）。"""
    from urllib.parse import urlparse, unquote

    parsed = urlparse(url)
    path = parsed.path.rstrip("/")
    name = os.path.basename(path)
    if not name:
        # 若无法解析出文件名，使用域名与时间戳占位
        host = parsed.netloc.replace(":", "_")
        name = f"download_{host}_{int(time.time())}"
    return unquote(name)


def download_with_requests(url: str, out_dir: Path, timeout: int, max_retry: int = 3) -> Tuple[bool, str, int]:
    """使用 requests 流式下载文件，返回 (成功标志, 保存路径, 字节数)。"""
    logging.info(f"准备下载文件: {url}")
    ensure_dir(out_dir)
    filename = infer_filename_from_url(url)
    out_path = out_dir / filename

    # 已存在且非空则跳过
    if out_path.exists() and out_path.stat().st_size > 0:
        logging.info("已存在且非空，跳过：%s", out_path)
        return True, str(out_path), out_path.stat().st_size

    for attempt in range(1, max_retry + 1):
        try:
            with requests.get(url, headers=DEFAULT_HEADERS, stream=True, timeout=timeout) as r:
                if r.status_code != 200:
                    raise RuntimeError(f"HTTP {r.status_code}")
                total = int(r.headers.get("Content-Length", 0))
                chunk = 1024 * 1024
                with open(out_path, "wb") as f, tqdm(
                    total=total if total > 0 else None,
                    unit="B",
                    unit_scale=True,
                    desc=f"Downloading {filename}",
                ) as pbar:
                    for part in r.iter_content(chunk_size=chunk):
                        if part:
                            f.write(part)
                            if total > 0:
                                pbar.update(len(part))
            size = out_path.stat().st_size
            logging.info("下载完成：%s | 大小：%s 字节", out_path, size)
            return True, str(out_path), size
        except Exception as e:  # noqa: E722
            logging.warning("第 %s 次下载失败：%s | 错误：%s", attempt, url, e)
            time.sleep(2 * attempt)

    logging.error("多次重试后仍失败：%s", url)
    return False, str(out_path), 0


def download_with_gdown(url: str, out_dir: Path, timeout: int) -> Tuple[bool, str, int]:
    """使用 gdown 下载 Google Drive 链接。"""
    logging.info(f"准备下载 Google Drive 文件: {url}")
    ensure_dir(out_dir)
    if gdown is None:
        logging.error("未安装 gdown，无法下载 Google Drive 链接：%s", url)
        return False, "", 0

    # 让 gdown 自动推断文件名；若失败则回退到占位名
    try:
        # gdown.download 支持输出路径，若为目录则会在目录下保存文件
        # 这里我们先切换工作目录以便保留原始文件名
        cwd = os.getcwd()
        os.chdir(str(out_dir))
        try:
            out = gdown.download(url=url, quiet=False)
        finally:
            os.chdir(cwd)
        if out is None:
            logging.error("gdown 返回空路径，下载失败：%s", url)
            return False, "", 0
        out_path = Path(out)
        size = out_path.stat().st_size if out_path.exists() else 0
        logging.info("下载完成(GDrive)：%s | 大小：%s 字节", out_path, size)
        return True, str(out_path), size
    except Exception as e:  # noqa: E722
        logging.error("gdown 下载失败：%s | 错误：%s", url, e)
        return False, "", 0


# ----------------------------- 主流程 -----------------------------

def write_manifest(manifest_path: Path, rows: List[Dict[str, str]]) -> None:
    """将下载清单写入 CSV 文件，便于后续溯源与制图。"""
    ensure_dir(manifest_path.parent)
    fieldnames = [
        "timestamp",
        "url",
        "category",
        "status",
        "saved_path",
        "bytes",
        "note",
    ]
    new_file = not manifest_path.exists()
    with open(manifest_path, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if new_file:
            writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    """主执行函数（参数已内置）。"""
    # -------- 参数设置 --------
    # 输出目录（绝对路径）
    out_dir = Path("E:/SDM01/data-gdw").resolve()
    # GDW 数据库主页 - 更新为更直接的 Figshare 数据仓库地址
    index_url = "https://figshare.com/articles/dataset/GDW_v1_0/25988293"
    # 抓取网页深度 (对于Figshare, 1层深度足够)
    max_depth = 1
    # 网络超时（秒）
    timeout = 60

    # -------- 路径准备 --------
    raw_dir = out_dir / "raw"
    log_dir = out_dir / "logs"
    manifest_path = out_dir / "manifest_gdw.csv"

    ensure_dir(out_dir)
    ensure_dir(raw_dir)
    setup_logging(log_dir)

    logging.info("开始抓取 GDW 链接 | index=%s | depth=%s", index_url, max_depth)
    direct_files, gdrive_files, visited_pages = crawl_and_collect(
        seed_urls=[index_url],
        max_depth=max_depth,
        timeout=timeout,
    )

    logging.info("抓取完成：%s 个网页 | %s 个直接文件 | %s 个GDrive", len(visited_pages), len(direct_files), len(gdrive_files))

    rows: List[Dict[str, str]] = []
    ts = dt.datetime.now().isoformat(timespec="seconds")

    # 记录已访问页面（便于审计）
    for pg in visited_pages:
        rows.append({
            "timestamp": ts,
            "url": pg,
            "category": "page",
            "status": "visited",
            "saved_path": "",
            "bytes": "",
            "note": "",
        })

    # 下载直接文件
    for url in direct_files:
        ok, save_path, nbytes = download_with_requests(url, raw_dir, timeout=timeout)
        rows.append({
            "timestamp": ts,
            "url": url,
            "category": "file",
            "status": "ok" if ok else "fail",
            "saved_path": save_path,
            "bytes": str(nbytes),
            "note": "requests",
        })

    # 下载 Google Drive 链接
    for url in gdrive_files:
        ok, save_path, nbytes = download_with_gdown(url, raw_dir, timeout=timeout)
        rows.append({
            "timestamp": ts,
            "url": url,
            "category": "gdrive",
            "status": "ok" if ok else "fail",
            "saved_path": save_path,
            "bytes": str(nbytes),
            "note": "gdown",
        })

    write_manifest(manifest_path, rows)
    logging.info("全部处理完成 | 清单：%s", manifest_path)


if __name__ == "__main__":
    main()


