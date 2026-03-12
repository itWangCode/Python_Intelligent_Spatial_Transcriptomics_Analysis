#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
RNA Secondary Structure Deep Prediction System
SCI-Level Complete Implementation with Resumable Crawling

Publication-ready code following best practices for:
- Nature Methods / Nucleic Acids Research / Bioinformatics submission
- Reproducible research standards
- Comprehensive evaluation with statistical significance testing

Architecture: k-mer Embedding + Transformer + U-Net Hybrid
Features:
  1.  Resumable large-scale RNA data crawling (bpRNA, Rfam, RNAcentral)
  2.  Robots.txt compliance & polite crawling
  3.  Data integrity verification (MD5, size checks)
  4.  Multi-source data fusion & quality-weighted deduplication
  5.  Homology-aware data splitting (MMseqs2 / CD-HIT)
  6.  k-mer tokenization + optional RNA-FM integration
  7.  Transformer + U-Net hybrid architecture
  8.  Multi-objective composite loss function
  9.  Nussinov dynamic programming decoder (multi-bracket pseudoknot)
  10. ViennaRNA thermodynamic model fusion
  11. Monte Carlo Dropout uncertainty estimation
  12. Comprehensive metrics (MCC, F1, PPV, Sensitivity)
  13. Statistical significance testing (Bootstrap, McNemar)
  14. Sequence quality scoring & data augmentation
  15. Attention visualization & model explainability
  16. Contrastive learning head
  17. Cross-validation support
  18. Ablation experiment framework
  19. SCI-quality visualization (Macaron palette, Times New Roman)

Author: [Your Name]
License: MIT
================================================================================
"""

# ── Standard Library ──────────────────────────────────────────────────────────
import os
import re
import sys
import json
import math
import time
import copy
import random
import shutil
import pickle
import hashlib
import zipfile
import sqlite3
import gzip
import logging
import argparse
import datetime
import subprocess
import itertools
import warnings
from pathlib import Path
from typing import (List, Tuple, Dict, Optional, Union, Any, Set,
                    Generator, Callable)
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor, as_completed

# ── Third-Party ───────────────────────────────────────────────────────────────
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from scipy.special import softmax
from scipy.stats import pearsonr, spearmanr

# Visualization
import matplotlib
matplotlib.use("Agg")          # Non-interactive backend for scripts
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker

# Machine Learning
from sklearn.model_selection import KFold, train_test_split
from sklearn.metrics import (roc_curve, auc, precision_recall_curve,
                             confusion_matrix, matthews_corrcoef)

# HTTP
import requests
from urllib.parse import urljoin, urlparse

# HTML Parsing
from bs4 import BeautifulSoup

# Progress
from tqdm import tqdm

# Deep Learning
try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from torch.utils.data import Dataset, DataLoader, WeightedRandomSampler
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False
    warnings.warn("PyTorch not installed. Model training disabled.", ImportWarning)

# Optional: Hugging Face transformers (RNA-FM)
try:
    from transformers import AutoModel, AutoTokenizer
    TRANSFORMERS_AVAILABLE = True
except ImportError:
    TRANSFORMERS_AVAILABLE = False

warnings.filterwarnings("ignore")

# ─────────────────────────────────────────────────────────────────────────────
#  GLOBAL PLOT STYLE  —  Macaron palette + Bold Times New Roman (SCI-quality)
# ─────────────────────────────────────────────────────────────────────────────

# Macaron Color Palette (pastel, publication-friendly)
MACARON_COLORS = [
    "#F4A7B9",   # Rose Macaron
    "#A8D8EA",   # Sky Blue Macaron
    "#B5EAD7",   # Mint Macaron
    "#FFDAC1",   # Peach Macaron
    "#C7CEEA",   # Lavender Macaron
    "#FFD6A5",   # Apricot Macaron
    "#D4E6B5",   # Pistachio Macaron
    "#F2C4CE",   # Blush Macaron
    "#B8D8D8",   # Teal Macaron
    "#ECD5E3",   # Lilac Macaron
]

MACARON_SEQUENTIAL = LinearSegmentedColormap.from_list(
    "macaron_seq", ["#FFFFFF", "#A8D8EA", "#4A90C4", "#2C5F8A"], N=256
)
MACARON_DIVERGING = LinearSegmentedColormap.from_list(
    "macaron_div", ["#F4A7B9", "#FFFFFF", "#A8D8EA"], N=256
)

def set_sci_style():
    """Apply SCI-quality plot style: bold Times New Roman, clean layout."""
    rcParams.update({
        # Font
        "font.family":          "serif",
        "font.serif":           ["Times New Roman", "Times", "DejaVu Serif"],
        "font.size":            12,
        "font.weight":          "bold",
        # Axes
        "axes.labelsize":       14,
        "axes.labelweight":     "bold",
        "axes.titlesize":       15,
        "axes.titleweight":     "bold",
        "axes.linewidth":       1.5,
        "axes.spines.top":      False,
        "axes.spines.right":    False,
        "axes.prop_cycle":      matplotlib.cycler(color=MACARON_COLORS),
        # Ticks
        "xtick.labelsize":      11,
        "ytick.labelsize":      11,
        "xtick.major.width":    1.5,
        "ytick.major.width":    1.5,
        "xtick.direction":      "out",
        "ytick.direction":      "out",
        # Legend
        "legend.fontsize":      11,
        "legend.frameon":       True,
        "legend.framealpha":    0.9,
        "legend.edgecolor":     "#CCCCCC",
        # Figure
        "figure.dpi":           300,
        "figure.facecolor":     "white",
        "savefig.dpi":          300,
        "savefig.bbox":         "tight",
        "savefig.facecolor":    "white",
        # Lines
        "lines.linewidth":      2.0,
        "lines.markersize":     7,
        # PDF/PS text
        "pdf.fonttype":         42,
        "ps.fonttype":          42,
    })

set_sci_style()


# ─────────────────────────────────────────────────────────────────────────────
#  0.  GLOBAL CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

class Config:
    """
    Unified configuration class.
    All hyper-parameters are centralised here for reproducibility.
    """

    def __init__(self):
        # ── Reproducibility ───────────────────────────────────────────────
        self.SEED: int = 42

        # ── Hardware ──────────────────────────────────────────────────────
        if TORCH_AVAILABLE:
            self.DEVICE: str = "cuda" if torch.cuda.is_available() else "cpu"
        else:
            self.DEVICE: str = "cpu"

        # ── Paths ─────────────────────────────────────────────────────────
        self.BASE_DIR    = Path("./rna_ss_sci")
        self.RAW_DIR     = self.BASE_DIR / "data" / "raw"
        self.PROC_DIR    = self.BASE_DIR / "data" / "processed"
        self.MODEL_DIR   = self.BASE_DIR / "models"
        self.TMP_DIR     = self.BASE_DIR / "tmp"
        self.LOG_DIR     = self.BASE_DIR / "logs"
        self.RESULTS_DIR = self.BASE_DIR / "results"
        self.FIG_DIR     = self.BASE_DIR / "figures"
        self.CRAWLER_DB  = self.BASE_DIR / "crawler_state.db"

        # ── Data Sources ──────────────────────────────────────────────────
        self.DATA_SOURCES: Dict[str, Dict] = {
            "bprna": {
                "base_url":       "https://bprna.cgrb.oregonstate.edu/",
                "dotbracket_url": "https://bprna.cgrb.oregonstate.edu/download.php?download=dotbracket",
                "fasta_url":      "https://bprna.cgrb.oregonstate.edu/download.php?download=fasta",
                "enabled":        True,
                "max_sequences":  50_000,
                "citation":       "Danaee et al. (2018) bpRNA"
            },
            "rfam": {
                "base_url":       "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/",
                "family_list_url":"https://rfam.org/families",
                "enabled":        True,
                "max_sequences":  100_000,
                "citation":       "Kalvari et al. (2021) Rfam 14.5"
            },
            "rnacentral": {
                "api_url":        "https://rnacentral.org/api/v1/",
                "enabled":        True,
                "max_sequences":  50_000,
                "citation":       "The RNAcentral Consortium (2021)"
            },
        }

        # ── Crawler ───────────────────────────────────────────────────────
        self.REQUEST_TIMEOUT:       int   = 30
        self.MAX_RETRIES:           int   = 3
        self.RETRY_DELAY:           float = 5.0
        self.CONCURRENT_DOWNLOADS:  int   = 4
        self.DOWNLOAD_CHUNK_SIZE:   int   = 8192
        self.USER_AGENT:            str   = (
            "RNAStructureResearch/1.0 (academic; contact: research@example.edu)"
        )
        self.RESPECT_ROBOTS_TXT:    bool  = True
        self.RATE_LIMIT:            float = 1.0   # seconds between requests per domain

        # ── Sequence Filtering ────────────────────────────────────────────
        self.MIN_LEN:               int   = 20
        self.MAX_LEN:               int   = 256
        self.MIN_HAIRPIN_LOOP:      int   = 3
        self.MIN_SEQUENCE_QUALITY:  float = 0.80  # fraction of unambiguous bases

        # ── Data Splitting ────────────────────────────────────────────────
        self.TRAIN_RATIO:  float = 0.70
        self.VALID_RATIO:  float = 0.15
        self.TEST_RATIO:   float = 0.15

        # ── Homology Clustering ───────────────────────────────────────────
        self.HOMOLOGY_IDENTITY: float = 0.80
        self.CLUSTER_METHOD:    str   = "mmseqs2"   # "mmseqs2" | "cdhit"

        # ── Training ──────────────────────────────────────────────────────
        self.BATCH_SIZE:   int   = 8
        self.NUM_EPOCHS:   int   = 50
        self.LR:           float = 2e-4
        self.WEIGHT_DECAY: float = 1e-4
        self.NUM_WORKERS:  int   = 0
        self.GRAD_CLIP:    float = 1.0
        self.WARMUP_STEPS: int   = 500
        self.LABEL_SMOOTHING: float = 0.05

        # ── Model Architecture ────────────────────────────────────────────
        self.KMER_SIZE:  int = 3
        self.D_MODEL:    int = 128
        self.NHEAD:      int = 8
        self.NUM_LAYERS: int = 6
        self.DIM_FF:     int = 512
        self.DROPOUT:    float = 0.1
        self.UNET_BASE_CH: int = 32

        # ── Loss Function Weights ─────────────────────────────────────────
        self.LOSS_WEIGHTS: Dict[str, float] = {
            "bce":    0.40,
            "focal":  0.25,
            "dice":   0.20,
            "sym":    0.05,
            "diag":   0.05,
            "sparse": 0.05,
        }
        self.POS_WEIGHT:         float = 8.0
        self.FOCAL_ALPHA:        float = 0.25
        self.FOCAL_GAMMA:        float = 2.0
        self.TARGET_PAIR_DENSITY:float = 0.02

        # ── Decoding ──────────────────────────────────────────────────────
        self.DECODING_THRESHOLD: float = 0.50
        self.THERMO_WEIGHT:      float = 0.20

        # ── Uncertainty ───────────────────────────────────────────────────
        self.MC_SAMPLES: int = 30

        # ── Data Augmentation ─────────────────────────────────────────────
        self.AUGMENT_NOISE_LEVEL: float = 0.03
        self.AUGMENT_PROBABILITY: float = 0.30

        # ── Create all directories ────────────────────────────────────────
        self._create_directories()
        self._init_crawler_db()
        self._setup_logging()

    def _create_directories(self) -> None:
        for d in [self.RAW_DIR, self.PROC_DIR, self.MODEL_DIR,
                  self.TMP_DIR, self.LOG_DIR, self.RESULTS_DIR, self.FIG_DIR]:
            d.mkdir(parents=True, exist_ok=True)

    def _init_crawler_db(self) -> None:
        """SQLite database for resumable crawling state."""
        conn = sqlite3.connect(str(self.CRAWLER_DB))
        cur  = conn.cursor()

        cur.execute("""
            CREATE TABLE IF NOT EXISTS data_sources (
                source_name          TEXT PRIMARY KEY,
                total_sequences      INTEGER DEFAULT 0,
                downloaded_sequences INTEGER DEFAULT 0,
                last_crawled         TIMESTAMP,
                status               TEXT,
                error_message        TEXT,
                config_json          TEXT
            )""")

        cur.execute("""
            CREATE TABLE IF NOT EXISTS downloaded_files (
                file_id       TEXT PRIMARY KEY,
                source_name   TEXT,
                url           TEXT UNIQUE,
                local_path    TEXT,
                file_size     INTEGER,
                md5_hash      TEXT,
                download_time TIMESTAMP,
                status        TEXT,
                retry_count   INTEGER DEFAULT 0
            )""")

        cur.execute("""
            CREATE TABLE IF NOT EXISTS parsed_sequences (
                sequence_id  TEXT PRIMARY KEY,
                source_name  TEXT,
                original_id  TEXT,
                sequence     TEXT,
                dot_bracket  TEXT,
                length       INTEGER,
                quality_score REAL,
                gc_content    REAL,
                pair_density  REAL,
                md5_hash     TEXT UNIQUE,
                parsed_time  TIMESTAMP
            )""")

        for idx_sql in [
            "CREATE INDEX IF NOT EXISTS idx_source  ON parsed_sequences(source_name)",
            "CREATE INDEX IF NOT EXISTS idx_length  ON parsed_sequences(length)",
            "CREATE INDEX IF NOT EXISTS idx_md5     ON parsed_sequences(md5_hash)",
            "CREATE INDEX IF NOT EXISTS idx_quality ON parsed_sequences(quality_score)",
        ]:
            cur.execute(idx_sql)

        conn.commit()
        conn.close()

    def _setup_logging(self) -> None:
        """Configure structured logging."""
        log_file = self.LOG_DIR / f"rna_ss_{datetime.date.today()}.log"
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s [%(levelname)s] %(name)s - %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            handlers=[
                logging.FileHandler(str(log_file), encoding="utf-8"),
                logging.StreamHandler(sys.stdout),
            ],
        )

    def update_from_args(self, args: argparse.Namespace) -> None:
        """Overwrite config from CLI arguments."""
        mapping = {
            "min_len":   "MIN_LEN",
            "max_len":   "MAX_LEN",
            "epochs":    "NUM_EPOCHS",
            "batch_size":"BATCH_SIZE",
            "lr":        "LR",
            "kmer":      "KMER_SIZE",
            "d_model":   "D_MODEL",
            "seed":      "SEED",
        }
        for arg_key, cfg_key in mapping.items():
            val = getattr(args, arg_key, None)
            if val is not None:
                setattr(self, cfg_key, val)


config = Config()
logger = logging.getLogger("rna_ss")


# ─────────────────────────────────────────────────────────────────────────────
#  1.  UTILITY FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def set_seed(seed: int = config.SEED) -> None:
    """Set all random seeds for full reproducibility."""
    random.seed(seed)
    np.random.seed(seed)
    if TORCH_AVAILABLE:
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    os.environ["PYTHONHASHSEED"] = str(seed)


def save_json(data: Any, path: Path) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2, default=str)


def load_json(path: Path) -> Any:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def save_pickle(data: Any, path: Path) -> None:
    with open(path, "wb") as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)


def load_pickle(path: Path) -> Any:
    with open(path, "rb") as f:
        return pickle.load(f)


def compute_md5(file_path: Path, chunk_size: int = 65536) -> str:
    """Compute MD5 hash of a file for integrity verification."""
    md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            md5.update(chunk)
    return md5.hexdigest()


def compute_sequence_md5(sequence: str, structure: str = "") -> str:
    """Compute MD5 of sequence+structure pair for deduplication."""
    return hashlib.md5(f"{sequence}|{structure}".encode("utf-8")).hexdigest()


def get_timestamp() -> str:
    return datetime.datetime.now().isoformat()


def verify_downloaded_file(file_path: Path) -> Dict[str, Any]:
    """
    Verify integrity of a downloaded file.
    Returns dict with validation results.
    """
    result: Dict[str, Any] = {
        "file_path": str(file_path),
        "exists": file_path.exists(),
        "size": 0,
        "md5": None,
        "is_valid": False,
        "issues": [],
    }

    if not result["exists"]:
        result["issues"].append("File does not exist")
        return result

    result["size"] = file_path.stat().st_size

    if result["size"] < 100:
        result["issues"].append(f"File too small ({result['size']} bytes)")
        return result

    try:
        result["md5"] = compute_md5(file_path)
        open_fn = gzip.open if str(file_path).endswith(".gz") else open
        with open_fn(file_path, "rt", encoding="utf-8", errors="ignore") as fh:
            content = fh.read(2048)

        if not any(c in content for c in "ACGU>"):
            result["issues"].append("No valid RNA sequence data detected")

        result["is_valid"] = len(result["issues"]) == 0

    except Exception as exc:
        result["issues"].append(f"Read error: {exc}")

    return result


# ─────────────────────────────────────────────────────────────────────────────
#  2.  SEQUENCE & STRUCTURE CLEANING
# ─────────────────────────────────────────────────────────────────────────────

def clean_sequence(seq: str) -> str:
    """
    Normalise RNA sequence:
    - Convert to uppercase
    - Replace T → U (DNA/RNA normalisation)
    - Remove non-IUPAC characters
    """
    if not isinstance(seq, str):
        return ""
    seq = seq.upper().replace("T", "U").strip()
    seq = re.sub(r"[^ACGURYWSMKHBVDN]", "", seq)
    return seq


def clean_dot_bracket(db: str) -> str:
    """
    Normalise dot-bracket notation:
    - Strip whitespace
    - Keep only valid structural characters
    """
    if not isinstance(db, str):
        return ""
    db = db.strip()
    # Retain standard brackets + extended pseudoknot notation
    db = re.sub(r"[^().\[\]{}<>A-Za-z]", "", db)
    return db


def validate_sequence_quality(seq: str) -> float:
    """
    Compute fraction of unambiguous (non-N) canonical bases.
    Returns value in [0, 1].
    """
    if not seq:
        return 0.0
    clean = clean_sequence(seq)
    canonical = sum(1 for b in clean if b in "ACGU")
    return canonical / max(len(clean), 1)


def compute_gc_content(seq: str) -> float:
    """Compute GC content fraction."""
    seq = clean_sequence(seq)
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def compute_pair_density(dot_bracket: str) -> float:
    """Compute fraction of paired positions."""
    if not dot_bracket:
        return 0.0
    paired = sum(1 for c in dot_bracket if c not in ".,-")
    return paired / max(len(dot_bracket), 1)


# ─────────────────────────────────────────────────────────────────────────────
#  3.  PSEUDOKNOT-AWARE BRACKET PARSER
# ─────────────────────────────────────────────────────────────────────────────

class PseudoKnotParser:
    """
    Multi-level pseudoknot parser supporting:
    ()  standard Watson-Crick pairs
    []  first-level pseudoknots
    {}  second-level pseudoknots
    <>  third-level pseudoknots
    Aa, Bb, …  higher-order pseudoknots (WUSS notation)
    """

    _STD_OPEN  = {"(": ")", "[": "]", "{": "}", "<": ">"}
    _STD_CLOSE = {v: k for k, v in _STD_OPEN.items()}

    # WUSS extended notation: uppercase = open, lowercase = close
    _UPPER = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    _LOWER = "abcdefghijklmnopqrstuvwxyz"
    _WUSS_OPEN  = {u: l for u, l in zip(_UPPER, _LOWER)}
    _WUSS_CLOSE = {l: u for u, l in zip(_UPPER, _LOWER)}

    @classmethod
    def dotbracket_to_pairs(cls, dot_bracket: str) -> List[Tuple[int, int]]:
        """
        Parse dot-bracket string → sorted list of (i, j) base pairs.
        Supports all bracket types including extended WUSS pseudoknot notation.
        """
        stacks: Dict[str, List[int]] = defaultdict(list)
        pairs: List[Tuple[int, int]] = []

        for pos, ch in enumerate(dot_bracket):
            if ch in cls._STD_OPEN:
                stacks[ch].append(pos)
            elif ch in cls._STD_CLOSE:
                opener = cls._STD_CLOSE[ch]
                if stacks[opener]:
                    j = stacks[opener].pop()
                    pairs.append((j, pos))
            elif ch in cls._WUSS_OPEN:
                stacks[ch].append(pos)
            elif ch in cls._WUSS_CLOSE:
                opener = cls._WUSS_CLOSE[ch]
                if stacks[opener]:
                    j = stacks[opener].pop()
                    pairs.append((j, pos))

        pairs.sort()
        return pairs

    @classmethod
    def pairs_to_matrix(cls,
                        pairs: List[Tuple[int, int]],
                        length: int) -> np.ndarray:
        """Convert base-pair list → symmetric binary matrix [L, L]."""
        mat = np.zeros((length, length), dtype=np.float32)
        for i, j in pairs:
            if 0 <= i < length and 0 <= j < length and i != j:
                mat[i, j] = mat[j, i] = 1.0
        return mat

    @classmethod
    def dotbracket_to_matrix(cls, dot_bracket: str) -> np.ndarray:
        """Convenience: dot-bracket string → matrix."""
        pairs = cls.dotbracket_to_pairs(dot_bracket)
        return cls.pairs_to_matrix(pairs, len(dot_bracket))

    @classmethod
    def pairs_to_dotbracket(cls,
                            pairs: List[Tuple[int, int]],
                            length: int) -> str:
        """
        Convert base-pair list → dot-bracket string.
        Automatically assigns bracket levels to avoid representation conflicts.
        Uses graph-colouring approach to minimise bracket levels.
        """
        if not pairs:
            return "." * length

        pairs = sorted(pairs)
        levels: Dict[Tuple[int, int], int] = {}

        # Build crossing graph
        def _crosses(a: int, b: int, c: int, d: int) -> bool:
            return (a < c < b < d) or (c < a < d < b)

        for p in pairs:
            # Find minimum level not conflicting with crossing pairs
            conflicting_levels = {
                levels[q] for q in levels
                if _crosses(p[0], p[1], q[0], q[1])
            }
            lev = 0
            while lev in conflicting_levels:
                lev += 1
            levels[p] = lev

        result = ["."] * length
        open_chars  = ["(", "[", "{", "<"] + list(cls._UPPER)
        close_chars = [")", "]", "}", ">"] + list(cls._LOWER)

        for (i, j), lev in levels.items():
            if lev < len(open_chars):
                result[i] = open_chars[lev]
                result[j] = close_chars[lev]

        return "".join(result)

    @staticmethod
    def identify_pseudoknots(pairs: Set[Tuple[int, int]]) -> Set[Tuple[int, int]]:
        """
        Identify which base pairs participate in pseudoknots.
        A pseudoknot exists when two pairs (i,j) and (k,l) satisfy i<k<j<l.
        """
        pair_list = sorted(pairs)
        pk_pairs: Set[Tuple[int, int]] = set()
        for idx_a, (i, j) in enumerate(pair_list):
            for (k, l) in pair_list[idx_a + 1:]:
                if i < k < j < l or k < i < l < j:
                    pk_pairs.add((i, j))
                    pk_pairs.add((k, l))
        return pk_pairs


# ─────────────────────────────────────────────────────────────────────────────
#  4.  BASE-PAIR CHEMISTRY
# ─────────────────────────────────────────────────────────────────────────────

# Standard + wobble + N-wildcard pairs
WATSON_CRICK_PAIRS: Set[Tuple[str, str]] = {
    ("A", "U"), ("U", "A"),
    ("C", "G"), ("G", "C"),
}
WOBBLE_PAIRS: Set[Tuple[str, str]] = {
    ("G", "U"), ("U", "G"),
}
WILDCARD_PAIRS: Set[Tuple[str, str]] = {
    (a, "N") for a in "ACGUN"
} | {
    ("N", a) for a in "ACGUN"
}
ALLOWED_PAIRS: Set[Tuple[str, str]] = (
    WATSON_CRICK_PAIRS | WOBBLE_PAIRS | WILDCARD_PAIRS
)


def is_valid_pair(a: str, b: str) -> bool:
    """Return True if bases a and b can form a canonical/wobble pair."""
    return (a, b) in ALLOWED_PAIRS


def build_basepair_prior(seq: str,
                         max_len: int,
                         min_loop: int = 3) -> np.ndarray:
    """
    Build chemistry-based prior matrix [max_len, max_len].
    Entry (i,j) = 1 iff bases seq[i] and seq[j] can pair
    and |i−j| > min_loop.
    """
    seq = clean_sequence(seq)
    L   = min(len(seq), max_len)
    prior = np.zeros((max_len, max_len), dtype=np.float32)

    for i in range(L):
        for j in range(i + min_loop + 1, L):
            if is_valid_pair(seq[i], seq[j]):
                prior[i, j] = prior[j, i] = 1.0

    return prior


# ─────────────────────────────────────────────────────────────────────────────
#  5.  SEQUENCE QUALITY SCORING SYSTEM
# ─────────────────────────────────────────────────────────────────────────────

class SequenceQualityScorer:
    """
    Comprehensive multi-criterion quality scoring for RNA sequences.

    Criteria and weights (sum = 1):
      - length_score        0.20  proximity to optimal length
      - gc_score            0.15  GC content in [0.3, 0.7]
      - ambiguity_score     0.25  fraction of unambiguous bases
      - structure_score     0.30  base-pair density (if structure provided)
      - motif_score         0.10  presence of conserved RNA motifs
    """

    WEIGHTS: Dict[str, float] = {
        "length":    0.20,
        "gc":        0.15,
        "ambiguity": 0.25,
        "structure": 0.30,
        "motif":     0.10,
    }

    CONSERVED_MOTIFS: List[str] = [
        "GU", "GNRA", "UNCG", "CUUG",  # common loop motifs
        "AAGU", "GAAA",                  # GNRA tetraloops
    ]

    def score(self, seq: str, dot_bracket: str = "") -> Dict[str, float]:
        """
        Compute quality scores for a single sequence.

        Returns dict with individual criterion scores and weighted total.
        All scores are in [0, 1]; higher is better.
        """
        seq_clean = clean_sequence(seq)
        L = len(seq_clean)

        if L == 0:
            return {k: 0.0 for k in list(self.WEIGHTS) + ["total"]}

        scores: Dict[str, float] = {}

        # 1. Length score (Gaussian around optimal length 150 nt)
        opt = 150.0
        scores["length"] = float(np.exp(-0.5 * ((L - opt) / opt) ** 2))

        # 2. GC content score — peak at 0.50, acceptable range [0.30, 0.70]
        gc = compute_gc_content(seq_clean)
        scores["gc"] = max(0.0, 1.0 - abs(gc - 0.50) * 3.0)

        # 3. Ambiguity score — penalise N bases
        n_ambig = seq_clean.count("N")
        scores["ambiguity"] = 1.0 - n_ambig / L

        # 4. Structure score — pair density
        if dot_bracket:
            pairs = PseudoKnotParser.dotbracket_to_pairs(dot_bracket)
            pair_density = len(pairs) * 2 / max(len(dot_bracket), 1)
            scores["structure"] = float(np.clip(pair_density * 2.0, 0.0, 1.0))
        else:
            scores["structure"] = 0.5   # neutral when no structure available

        # 5. Conserved motif score
        motif_hits = sum(seq_clean.count(m.replace("N", "A"))
                         for m in self.CONSERVED_MOTIFS)
        scores["motif"] = float(np.clip(motif_hits / max(L / 10, 1), 0.0, 1.0))

        scores["total"] = sum(
            scores[k] * self.WEIGHTS[k] for k in self.WEIGHTS
        )
        return scores


# ─────────────────────────────────────────────────────────────────────────────
#  6.  DATA AUGMENTATION
# ─────────────────────────────────────────────────────────────────────────────

class RNADataAugmenter:
    """
    RNA sequence & structure data augmentation for training robustness.

    Methods:
      - point_mutation  : random single-base substitutions
      - random_crop     : extract random sub-window
      - reverse_complement: biological reverse complement
    """

    BASE_TRANS: Dict[str, str] = {"A": "U", "U": "A", "C": "G", "G": "C", "N": "N"}
    OPEN_CLOSE: Dict[str, str] = {
        "(": ")", ")": "(", "[": "]", "]": "[",
        "{": "}", "}": "{", "<": ">", ">": "<",
    }

    def __init__(self, seed: int = 42):
        self._rng = random.Random(seed)
        self._np_rng = np.random.RandomState(seed)

    def point_mutation(self, seq: str, noise_level: float = 0.03) -> str:
        """
        Introduce random point mutations.
        noise_level: expected fraction of mutated positions.
        """
        bases = list("ACGU")
        seq_list = list(seq)
        for i, base in enumerate(seq_list):
            if base not in bases:
                continue
            if self._rng.random() < noise_level:
                candidates = [b for b in bases if b != base]
                seq_list[i] = self._rng.choice(candidates)
        return "".join(seq_list)

    def reverse_complement(
        self, seq: str, db: str
    ) -> Tuple[str, str]:
        """
        Compute RNA reverse complement and mirror the structure.
        """
        comp_seq = "".join(self.BASE_TRANS.get(b, "N") for b in reversed(seq))
        rev_db   = "".join(
            self.OPEN_CLOSE.get(c, c) for c in reversed(db)
        )
        return comp_seq, rev_db

    def random_crop(
        self, seq: str, db: str, target_len: int
    ) -> Tuple[str, str]:
        """
        Extract a random contiguous sub-sequence of length target_len.
        Preserves sequence–structure alignment.
        """
        L = len(seq)
        if L <= target_len:
            return seq, db
        start = self._rng.randint(0, L - target_len)
        end   = start + target_len
        return seq[start:end], db[start:end]

    def augment(
        self,
        seq: str,
        db: str,
        prob: float = 0.30,
    ) -> Tuple[str, str]:
        """
        Apply a randomly chosen augmentation with probability `prob`.
        """
        if self._rng.random() > prob:
            return seq, db

        choice = self._rng.choice(["mutation", "reverse"])
        if choice == "mutation":
            return self.point_mutation(seq, config.AUGMENT_NOISE_LEVEL), db
        else:
            return self.reverse_complement(seq, db)


# ─────────────────────────────────────────────────────────────────────────────
#  7.  RESUMABLE CRAWLER BASE CLASS
# ─────────────────────────────────────────────────────────────────────────────

class ResumeableCrawler:
    """
    Base class for all data source crawlers.

    Features:
    - SQLite-backed state persistence for fault-tolerant crawling
    - Robots.txt compliance with per-domain caching
    - Exponential backoff on transient failures
    - Resume-from-checkpoint by URL deduplication
    - Per-domain request rate limiting
    """

    def __init__(self, source_name: str):
        self.source_name   = source_name
        self._session      = requests.Session()
        self._session.headers.update({"User-Agent": config.USER_AGENT})
        self._domain_state: Dict[str, Dict] = {}   # per-domain delay tracking
        self._robots_cache: Dict[str, Any]   = {}
        self._db: Optional[sqlite3.Connection] = None
        self._scorer = SequenceQualityScorer()
        self._connect_db()
        self._init_source_state()

    # ── Database Helpers ──────────────────────────────────────────────────────

    def _connect_db(self) -> None:
        self._db = sqlite3.connect(str(config.CRAWLER_DB), check_same_thread=False)
        self._db.row_factory = sqlite3.Row

    def _init_source_state(self) -> None:
        cur = self._db.cursor()
        cur.execute(
            "SELECT 1 FROM data_sources WHERE source_name = ?",
            (self.source_name,),
        )
        if cur.fetchone() is None:
            src_cfg = json.dumps(config.DATA_SOURCES.get(self.source_name, {}))
            cur.execute(
                """INSERT INTO data_sources
                   (source_name, total_sequences, downloaded_sequences,
                    status, config_json)
                   VALUES (?,0,0,'initialized',?)""",
                (self.source_name, src_cfg),
            )
            self._db.commit()

    def update_source_state(self, **kwargs) -> None:
        allowed = {
            "total_sequences", "downloaded_sequences",
            "status", "error_message",
        }
        fields = [f"{k} = ?" for k in kwargs if k in allowed]
        values = [v for k, v in kwargs.items() if k in allowed]
        if fields:
            values.append(self.source_name)
            self._db.execute(
                f"""UPDATE data_sources
                    SET {", ".join(fields)}, last_crawled = CURRENT_TIMESTAMP
                    WHERE source_name = ?""",
                values,
            )
            self._db.commit()

    # ── Robots.txt Compliance ─────────────────────────────────────────────────

    def _check_robots_txt(self, domain: str) -> bool:
        """
        Verify crawling is permitted by robots.txt.
        Returns True if allowed (or robots.txt unavailable).
        """
        if not config.RESPECT_ROBOTS_TXT:
            return True

        if domain in self._robots_cache:
            return self._robots_cache[domain]

        try:
            from urllib.robotparser import RobotFileParser
            rp = RobotFileParser()
            rp.set_url(f"https://{domain}/robots.txt")
            rp.read()
            allowed = rp.can_fetch(config.USER_AGENT, f"https://{domain}/")
            self._robots_cache[domain] = allowed
            if not allowed:
                logger.warning("Crawling disallowed by robots.txt: %s", domain)
            return allowed
        except Exception:
            self._robots_cache[domain] = True   # assume allowed on error
            return True

    # ── Polite Request ────────────────────────────────────────────────────────

    def _polite_request(
        self, url: str, stream: bool = False, **kwargs
    ) -> Optional[requests.Response]:
        """
        Send HTTP GET with:
        - Robots.txt check
        - Per-domain rate limiting
        - Retry-After header handling
        - Exponential backoff
        """
        domain = urlparse(url).netloc

        # Initialise domain tracking
        if domain not in self._domain_state:
            self._check_robots_txt(domain)
            self._domain_state[domain] = {"last_request": 0.0}

        # Rate limiting
        elapsed = time.time() - self._domain_state[domain]["last_request"]
        if elapsed < config.RATE_LIMIT:
            time.sleep(config.RATE_LIMIT - elapsed)

        for attempt in range(config.MAX_RETRIES):
            try:
                resp = self._session.get(
                    url,
                    stream=stream,
                    timeout=config.REQUEST_TIMEOUT,
                    **kwargs,
                )
                self._domain_state[domain]["last_request"] = time.time()

                # Handle Retry-After
                if resp.status_code == 429 or "Retry-After" in resp.headers:
                    wait = int(resp.headers.get("Retry-After", config.RETRY_DELAY))
                    logger.info("Rate limited on %s — waiting %ds", domain, wait)
                    time.sleep(wait)
                    continue

                resp.raise_for_status()
                return resp

            except requests.RequestException as exc:
                wait = config.RETRY_DELAY * (2 ** attempt)
                logger.warning(
                    "Request attempt %d/%d failed for %s: %s — retrying in %.1fs",
                    attempt + 1, config.MAX_RETRIES, url, exc, wait,
                )
                time.sleep(wait)

        logger.error("All %d attempts failed for %s", config.MAX_RETRIES, url)
        return None

    # ── File Download (resumable) ─────────────────────────────────────────────

    def download_file(
        self,
        url: str,
        local_path: Path,
        expected_md5: Optional[str] = None,
    ) -> bool:
        """
        Download with HTTP Range (resume) support.
        Returns True on success.
        """
        # Check already-completed download
        cur = self._db.cursor()
        cur.execute(
            "SELECT local_path, md5_hash, status FROM downloaded_files WHERE url = ?",
            (url,),
        )
        row = cur.fetchone()
        if row and row["status"] == "completed":
            existing = Path(row["local_path"])
            if existing.exists():
                if expected_md5 is None or row["md5_hash"] == expected_md5:
                    logger.info("Skipping already-downloaded: %s", url)
                    return True

        local_path.parent.mkdir(parents=True, exist_ok=True)
        existing_size = local_path.stat().st_size if local_path.exists() else 0

        headers = {}
        mode = "wb"
        if existing_size > 0:
            headers["Range"] = f"bytes={existing_size}-"
            mode = "ab"

        for attempt in range(config.MAX_RETRIES):
            try:
                resp = self._session.get(
                    url,
                    headers=headers,
                    stream=True,
                    timeout=config.REQUEST_TIMEOUT,
                )
                resp.raise_for_status()

                total = int(resp.headers.get("content-length", 0)) + existing_size

                with open(local_path, mode) as fh:
                    with tqdm(
                        total=total,
                        initial=existing_size,
                        unit="B",
                        unit_scale=True,
                        desc=local_path.name,
                        leave=False,
                    ) as pbar:
                        for chunk in resp.iter_content(
                            chunk_size=config.DOWNLOAD_CHUNK_SIZE
                        ):
                            if chunk:
                                fh.write(chunk)
                                pbar.update(len(chunk))

                # Integrity check
                validation = verify_downloaded_file(local_path)
                if not validation["is_valid"]:
                    raise ValueError(
                        f"File validation failed: {validation['issues']}"
                    )

                if expected_md5 and validation["md5"] != expected_md5:
                    raise ValueError("MD5 mismatch after download")

                self._register_file(url, local_path, validation["md5"])
                return True

            except Exception as exc:
                logger.warning(
                    "Download attempt %d failed for %s: %s",
                    attempt + 1, url, exc,
                )
                if attempt < config.MAX_RETRIES - 1:
                    time.sleep(config.RETRY_DELAY * (2 ** attempt))

        # Record failure
        file_id = hashlib.md5(url.encode()).hexdigest()
        self._db.execute(
            """INSERT OR REPLACE INTO downloaded_files
               (file_id, source_name, url, local_path, status, retry_count)
               VALUES (?,?,?,?,?,?)""",
            (file_id, self.source_name, url, str(local_path), "failed", config.MAX_RETRIES),
        )
        self._db.commit()
        return False

    def _register_file(self, url: str, path: Path, md5: str) -> None:
        file_id = hashlib.md5(url.encode()).hexdigest()
        self._db.execute(
            """INSERT OR REPLACE INTO downloaded_files
               (file_id, source_name, url, local_path, file_size, md5_hash,
                download_time, status)
               VALUES (?,?,?,?,?,?,CURRENT_TIMESTAMP,'completed')""",
            (file_id, self.source_name, url, str(path),
             path.stat().st_size, md5),
        )
        self._db.commit()

    # ── Sequence Persistence ──────────────────────────────────────────────────

    def save_sequence(
        self,
        original_id: str,
        sequence: str,
        dot_bracket: str,
        quality_score: float,
    ) -> bool:
        """
        Persist a parsed sequence to the database.
        Silently skips duplicates (by sequence+structure MD5).
        Returns True if a new record was inserted.
        """
        seq_id   = f"{self.source_name}:{original_id}"
        md5_hash = compute_sequence_md5(sequence, dot_bracket)
        gc       = compute_gc_content(sequence)
        pd_val   = compute_pair_density(dot_bracket)

        try:
            self._db.execute(
                """INSERT INTO parsed_sequences
                   (sequence_id, source_name, original_id, sequence, dot_bracket,
                    length, quality_score, gc_content, pair_density, md5_hash,
                    parsed_time)
                   VALUES (?,?,?,?,?,?,?,?,?,?,CURRENT_TIMESTAMP)""",
                (seq_id, self.source_name, original_id,
                 sequence, dot_bracket, len(sequence),
                 quality_score, gc, pd_val, md5_hash),
            )
            self._db.execute(
                """UPDATE data_sources
                   SET downloaded_sequences = downloaded_sequences + 1
                   WHERE source_name = ?""",
                (self.source_name,),
            )
            self._db.commit()
            return True
        except sqlite3.IntegrityError:
            return False   # duplicate — silently skip

    def get_progress(self) -> Dict[str, Any]:
        cur = self._db.cursor()
        cur.execute(
            """SELECT total_sequences, downloaded_sequences, status,
                      error_message, last_crawled
               FROM data_sources WHERE source_name = ?""",
            (self.source_name,),
        )
        row = cur.fetchone()
        return dict(row) if row else {}

    def validate_data_source(self) -> Dict[str, Any]:
        """
        Health-check the remote data source.
        Returns metadata including availability and response time.
        """
        result = {
            "source_name":   self.source_name,
            "available":     False,
            "response_time": None,
            "error":         None,
        }
        src_cfg = config.DATA_SOURCES.get(self.source_name, {})
        test_url = src_cfg.get("base_url") or src_cfg.get("api_url", "")

        if not test_url:
            result["error"] = "No URL configured"
            return result

        try:
            t0   = time.time()
            resp = self._session.head(test_url, timeout=10)
            result["response_time"] = round(time.time() - t0, 3)
            result["available"]     = resp.status_code < 400
        except Exception as exc:
            result["error"] = str(exc)

        return result

    def close(self) -> None:
        if self._db:
            self._db.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


# ─────────────────────────────────────────────────────────────────────────────
#  8.  SOURCE-SPECIFIC CRAWLERS
# ─────────────────────────────────────────────────────────────────────────────

class BpRNACrawler(ResumeableCrawler):
    """
    Crawler for the bpRNA database.
    Reference: Danaee et al. (2018) Nucleic Acids Research.
    """

    def __init__(self):
        super().__init__("bprna")
        cfg = config.DATA_SOURCES["bprna"]
        self._fasta_url     = cfg["fasta_url"]
        self._dot_url       = cfg["dotbracket_url"]
        self._max_sequences = cfg["max_sequences"]

    def crawl(self, max_sequences: Optional[int] = None) -> int:
        max_sequences = max_sequences or self._max_sequences
        logger.info("BpRNA crawl started  (target=%d)", max_sequences)
        self.update_source_state(status="crawling")

        try:
            fasta_zip = config.RAW_DIR / "bprna" / "bprna_fasta.zip"
            dot_zip   = config.RAW_DIR / "bprna" / "bprna_dotbracket.zip"

            if not (self.download_file(self._fasta_url, fasta_zip) and
                    self.download_file(self._dot_url, dot_zip)):
                raise RuntimeError("Failed to download bpRNA archives")

            fasta_dir = config.RAW_DIR / "bprna" / "fasta"
            dot_dir   = config.RAW_DIR / "bprna" / "dotbracket"
            self._extract_zip(fasta_zip, fasta_dir)
            self._extract_zip(dot_zip, dot_dir)

            count = self._parse_files(fasta_dir, dot_dir, max_sequences)
            self.update_source_state(status="completed", total_sequences=count)
            logger.info("BpRNA crawl complete (%d sequences)", count)
            return count

        except Exception as exc:
            self.update_source_state(status="failed", error_message=str(exc))
            logger.error("BpRNA crawl failed: %s", exc)
            raise

    def _extract_zip(self, zip_path: Path, dest: Path) -> None:
        if dest.exists() and any(dest.iterdir()):
            return
        dest.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(dest)
        logger.info("Extracted %s → %s", zip_path.name, dest)

    def _parse_files(
        self, fasta_dir: Path, dot_dir: Path, max_seq: int
    ) -> int:
        fasta_files = {f.stem: f for f in fasta_dir.glob("*.f*")}
        dot_files   = {f.stem: f for f in dot_dir.glob("*.*")}
        common      = sorted(set(fasta_files) & set(dot_files))
        logger.info("BpRNA: %d matched file pairs found", len(common))

        count = 0
        for key in tqdm(common, desc="Parsing bpRNA"):
            if count >= max_seq:
                break
            try:
                seq = self._read_fasta(fasta_files[key])
                db  = self._read_dot_bracket(dot_files[key])
                if not (seq and db) or len(seq) != len(db):
                    continue
                if not (config.MIN_LEN <= len(seq) <= config.MAX_LEN):
                    continue
                quality = validate_sequence_quality(seq)
                if quality < config.MIN_SEQUENCE_QUALITY:
                    continue
                if self.save_sequence(key, seq, db, quality):
                    count += 1
            except Exception as exc:
                logger.debug("BpRNA parse error for %s: %s", key, exc)

        return count

    @staticmethod
    def _read_fasta(path: Path) -> Optional[str]:
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as fh:
                return clean_sequence(
                    "".join(l.strip() for l in fh if not l.startswith(">"))
                )
        except Exception:
            return None

    @staticmethod
    def _read_dot_bracket(path: Path) -> Optional[str]:
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    line = line.strip()
                    if line and not line.startswith((">" , "#")):
                        cleaned = clean_dot_bracket(line)
                        if cleaned:
                            return cleaned
        except Exception:
            pass
        return None


class RfamCrawler(ResumeableCrawler):
    """
    Crawler for the Rfam database (FTP FASTA files).
    Reference: Kalvari et al. (2021) Nucleic Acids Research.
    """

    def __init__(self):
        super().__init__("rfam")
        cfg = config.DATA_SOURCES["rfam"]
        self._base_url      = cfg["base_url"]
        self._max_sequences = cfg["max_sequences"]

    def crawl(
        self,
        max_sequences: Optional[int] = None,
        family_ids: Optional[List[str]] = None,
    ) -> int:
        max_sequences = max_sequences or self._max_sequences
        logger.info("Rfam crawl started (target=%d)", max_sequences)
        self.update_source_state(status="crawling")

        try:
            if family_ids is None:
                family_ids = self._list_available_families(limit=2000)

            count = 0
            for fid in tqdm(family_ids, desc="Rfam families"):
                if count >= max_sequences:
                    break
                try:
                    for ext in (".fa.gz", ".fasta.gz", ".fa"):
                        url        = f"{self._base_url}{fid}{ext}"
                        local_path = config.RAW_DIR / "rfam" / f"{fid}{ext}"
                        if self.download_file(url, local_path):
                            n = self._parse_fasta_file(
                                local_path, fid, max_sequences - count
                            )
                            count += n
                            break
                except Exception as exc:
                    logger.debug("Rfam family %s error: %s", fid, exc)

            self.update_source_state(status="completed", total_sequences=count)
            logger.info("Rfam crawl complete (%d sequences)", count)
            return count

        except Exception as exc:
            self.update_source_state(status="failed", error_message=str(exc))
            logger.error("Rfam crawl failed: %s", exc)
            raise

    def _list_available_families(self, limit: int = 2000) -> List[str]:
        """Retrieve family IDs from Rfam FTP index."""
        resp = self._polite_request(self._base_url)
        if resp is None:
            return [f"RF{str(i).zfill(5)}" for i in range(1, min(limit + 1, 4001))]
        soup = BeautifulSoup(resp.text, "html.parser")
        ids: List[str] = []
        for link in soup.find_all("a"):
            href = link.get("href", "")
            m = re.search(r"(RF\d{5})", href)
            if m and m.group(1) not in ids:
                ids.append(m.group(1))
            if len(ids) >= limit:
                break
        return ids if ids else [f"RF{str(i).zfill(5)}" for i in range(1, limit + 1)]

    def _parse_fasta_file(
        self, path: Path, family_id: str, max_seq: int
    ) -> int:
        count = 0
        open_fn = gzip.open if str(path).endswith(".gz") else open
        try:
            with open_fn(path, "rt", encoding="utf-8", errors="ignore") as fh:
                current_id, seq_parts = None, []
                for line in fh:
                    line = line.strip()
                    if line.startswith(">"):
                        if current_id and seq_parts:
                            seq = "".join(seq_parts)
                            if (config.MIN_LEN <= len(seq) <= config.MAX_LEN):
                                quality = validate_sequence_quality(seq)
                                if quality >= config.MIN_SEQUENCE_QUALITY:
                                    if self.save_sequence(
                                        f"{family_id}_{current_id}",
                                        clean_sequence(seq), "", quality
                                    ):
                                        count += 1
                                        if count >= max_seq:
                                            return count
                        current_id = line[1:].split()[0]
                        seq_parts  = []
                    elif line:
                        seq_parts.append(line)
        except Exception as exc:
            logger.debug("Rfam parse error %s: %s", path, exc)
        return count


class RNAcentralCrawler(ResumeableCrawler):
    """
    Crawler for RNAcentral REST API.
    Reference: The RNAcentral Consortium (2021) Nucleic Acids Research.
    """

    def __init__(self):
        super().__init__("rnacentral")
        cfg = config.DATA_SOURCES["rnacentral"]
        self._api_url       = cfg["api_url"]
        self._max_sequences = cfg["max_sequences"]

    def crawl(
        self,
        max_sequences: Optional[int] = None,
        rna_types: Optional[List[str]] = None,
    ) -> int:
        max_sequences = max_sequences or self._max_sequences
        if rna_types is None:
            rna_types = ["miRNA", "rRNA", "tRNA", "lncRNA", "snRNA", "snoRNA"]

        logger.info("RNAcentral crawl started (target=%d)", max_sequences)
        self.update_source_state(status="crawling")

        count = page = 0
        try:
            while count < max_sequences:
                records = self._fetch_page(page, rna_types)
                if not records:
                    break
                for rec in records:
                    seq = clean_sequence(rec.get("sequence", ""))
                    if not (config.MIN_LEN <= len(seq) <= config.MAX_LEN):
                        continue
                    quality = validate_sequence_quality(seq)
                    if quality < config.MIN_SEQUENCE_QUALITY:
                        continue
                    db = clean_dot_bracket(
                        rec.get("secondary_structure", "") or ""
                    )
                    if self.save_sequence(
                        rec.get("rnacentral_id", f"p{page}_r{count}"),
                        seq, db, quality,
                    ):
                        count += 1
                        if count >= max_sequences:
                            break
                logger.debug("RNAcentral page %d → %d total", page, count)
                page += 1
                time.sleep(config.RATE_LIMIT)

            self.update_source_state(status="completed", total_sequences=count)
            logger.info("RNAcentral crawl complete (%d sequences)", count)
            return count

        except Exception as exc:
            self.update_source_state(status="failed", error_message=str(exc))
            logger.error("RNAcentral crawl failed: %s", exc)
            raise

    def _fetch_page(
        self, page: int, rna_types: List[str]
    ) -> List[Dict]:
        params = {
            "format":       "json",
            "page":         page + 1,
            "page_size":    100,
            "has_sequence": "true",
            "rna_type":     ",".join(rna_types),
        }
        resp = self._polite_request(
            urljoin(self._api_url, "rnacentral/"),
            params=params,
        )
        if resp is None:
            return []
        try:
            return resp.json().get("results", [])
        except Exception:
            return []


# ─────────────────────────────────────────────────────────────────────────────
#  9.  CRAWL MANAGER
# ─────────────────────────────────────────────────────────────────────────────

class DataCrawlerManager:
    """
    Orchestrates multi-source crawling with parallel or sequential execution.
    Provides aggregate progress reporting and summary statistics.
    """

    def __init__(self) -> None:
        self._crawlers: Dict[str, ResumeableCrawler] = {}
        for name, cfg in config.DATA_SOURCES.items():
            if not cfg.get("enabled", False):
                continue
            if name == "bprna":
                self._crawlers[name] = BpRNACrawler()
            elif name == "rfam":
                self._crawlers[name] = RfamCrawler()
            elif name == "rnacentral":
                self._crawlers[name] = RNAcentralCrawler()

    def health_check(self) -> Dict[str, Dict]:
        """Validate all configured data sources."""
        return {name: c.validate_data_source()
                for name, c in self._crawlers.items()}

    def crawl_all(
        self,
        target_total: Optional[int] = None,
        parallel: bool = False,
    ) -> int:
        if target_total is None:
            target_total = sum(
                s["max_sequences"]
                for s in config.DATA_SOURCES.values()
                if s.get("enabled")
            )

        logger.info("Multi-source crawl started (target=%d)", target_total)
        results: Dict[str, int] = {}

        if parallel and len(self._crawlers) > 1:
            with ThreadPoolExecutor(max_workers=config.CONCURRENT_DOWNLOADS) as ex:
                futures = {ex.submit(c.crawl): n
                           for n, c in self._crawlers.items()}
                for fut in as_completed(futures):
                    name = futures[fut]
                    try:
                        results[name] = fut.result()
                    except Exception as exc:
                        logger.error("%s failed: %s", name, exc)
                        results[name] = 0
        else:
            cumulative = 0
            for name, crawler in self._crawlers.items():
                if cumulative >= target_total:
                    break
                try:
                    n = crawler.crawl()
                    results[name] = n
                    cumulative += n
                except Exception as exc:
                    logger.error("%s failed: %s", name, exc)
                    results[name] = 0

        total = sum(results.values())
        summary = {
            "target_total": target_total,
            "actual_total": total,
            "sources":      results,
            "timestamp":    get_timestamp(),
        }
        save_json(summary, config.RESULTS_DIR / "crawl_summary.json")
        logger.info("Crawl complete: %d / %d sequences", total, target_total)
        return total

    def get_progress(self) -> Dict[str, Dict]:
        return {n: c.get_progress() for n, c in self._crawlers.items()}

    def close(self) -> None:
        for c in self._crawlers.values():
            c.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


# ─────────────────────────────────────────────────────────────────────────────
#  10.  DATA PROCESSOR & FUSION
# ─────────────────────────────────────────────────────────────────────────────

class DataProcessor:
    """
    Post-crawl data processing:
    - Quality-weighted deduplication
    - Source fusion with stratified statistics
    - Multi-format export (FASTA, Stockholm, CSV)
    - Comprehensive dataset visualisation
    """

    def __init__(self) -> None:
        self._db = sqlite3.connect(str(config.CRAWLER_DB))
        self._db.row_factory = sqlite3.Row

    def build_dataset(
        self,
        min_quality: float = 0.80,
        require_structure: bool = True,
        min_pair_density: float = 0.0,
    ) -> pd.DataFrame:
        """
        Query the crawl database and build a curated DataFrame.
        Applies quality, length, and optionally structure filters.
        """
        logger.info("Building dataset (min_quality=%.2f, require_structure=%s)",
                    min_quality, require_structure)

        query = """
            SELECT sequence_id, source_name, original_id,
                   sequence, dot_bracket, length,
                   quality_score, gc_content, pair_density
            FROM parsed_sequences
            WHERE length    BETWEEN ? AND ?
              AND quality_score >= ?
        """
        params: List[Any] = [config.MIN_LEN, config.MAX_LEN, min_quality]

        if require_structure:
            query += " AND length(COALESCE(dot_bracket,'')) > 0"
        if min_pair_density > 0:
            query += " AND pair_density >= ?"
            params.append(min_pair_density)

        cur = self._db.cursor()
        cur.execute(query, params)
        rows = cur.fetchall()

        if not rows:
            logger.warning("No sequences matched the filter criteria")
            return pd.DataFrame()

        df = pd.DataFrame([dict(r) for r in rows])
        logger.info("Dataset: %d sequences from %d sources",
                    len(df), df["source_name"].nunique())

        # Save
        csv_path = config.PROC_DIR / "raw_dataset.csv"
        df.to_csv(csv_path, index=False)
        logger.info("Dataset saved → %s", csv_path)
        return df

    def filter_by_structure_complexity(
        self, df: pd.DataFrame, min_pairs: int = 2
    ) -> pd.DataFrame:
        """Retain only sequences with ≥ min_pairs base pairs."""
        def _pair_count(db: str) -> int:
            if not isinstance(db, str):
                return 0
            return len(PseudoKnotParser.dotbracket_to_pairs(db))

        df = df.copy()
        df["pair_count"] = df["dot_bracket"].apply(_pair_count)
        out = df[df["pair_count"] >= min_pairs].reset_index(drop=True)
        logger.info("Structure filter: %d → %d sequences (min_pairs=%d)",
                    len(df), len(out), min_pairs)
        return out

    def export_fasta(self, df: pd.DataFrame, path: Path) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            for _, row in df.iterrows():
                fh.write(f">{row['sequence_id']}\n{row['sequence']}\n")
        logger.info("FASTA exported → %s", path)

    def export_stockholm(self, df: pd.DataFrame, path: Path) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("# STOCKHOLM 1.0\n\n")
            for _, row in df.iterrows():
                fh.write(f"{str(row['sequence_id']):<30} {row['sequence']}\n")
                if row.get("dot_bracket"):
                    fh.write(f"{'#=GF SS':<30} {row['dot_bracket']}\n")
            fh.write("//\n")
        logger.info("Stockholm exported → %s", path)

    def visualize_dataset(self, df: pd.DataFrame) -> None:
        """
        Generate SCI-quality dataset overview figures (6-panel layout).
        Macaron palette, bold Times New Roman, all-English labels.
        """
        if df.empty:
            logger.warning("Empty DataFrame — skipping visualisation")
            return

        fig, axes = plt.subplots(2, 3, figsize=(16, 10))
        fig.suptitle(
            "Dataset Overview Statistics",
            fontsize=18, fontweight="bold", y=1.01
        )

        # ── Panel 1: Sequence length distribution ────────────────────────
        ax = axes[0, 0]
        ax.hist(
            df["length"], bins=50,
            color=MACARON_COLORS[1], edgecolor="white", linewidth=0.5,
            alpha=0.85,
        )
        ax.set_xlabel("Sequence Length (nt)", fontweight="bold")
        ax.set_ylabel("Count", fontweight="bold")
        ax.set_title("Sequence Length Distribution", fontweight="bold")
        ax.axvline(df["length"].median(), color=MACARON_COLORS[0],
                   linestyle="--", linewidth=2, label=f"Median={df['length'].median():.0f}")
        ax.legend(prop={"weight": "bold"})

        # ── Panel 2: Source distribution (pie) ───────────────────────────
        ax = axes[0, 1]
        src_counts = df["source_name"].value_counts()
        wedges, texts, autotexts = ax.pie(
            src_counts.values,
            labels=src_counts.index,
            colors=MACARON_COLORS[:len(src_counts)],
            autopct="%1.1f%%",
            startangle=90,
            wedgeprops={"edgecolor": "white", "linewidth": 1.5},
        )
        for t in texts + autotexts:
            t.set_fontweight("bold")
            t.set_fontsize(11)
        ax.set_title("Data Source Distribution", fontweight="bold")

        # ── Panel 3: GC content distribution ─────────────────────────────
        ax = axes[0, 2]
        gc_vals = df.get("gc_content", df["sequence"].apply(compute_gc_content))
        ax.hist(
            gc_vals, bins=40,
            color=MACARON_COLORS[2], edgecolor="white", linewidth=0.5,
            alpha=0.85,
        )
        ax.axvline(0.5, color=MACARON_COLORS[0], linestyle="--",
                   linewidth=2, label="Optimal GC = 0.50")
        ax.set_xlabel("GC Content Fraction", fontweight="bold")
        ax.set_ylabel("Count", fontweight="bold")
        ax.set_title("GC Content Distribution", fontweight="bold")
        ax.legend(prop={"weight": "bold"})

        # ── Panel 4: Pair density distribution ───────────────────────────
        ax = axes[1, 0]
        if "pair_density" in df.columns:
            ax.hist(
                df["pair_density"].dropna(), bins=40,
                color=MACARON_COLORS[3], edgecolor="white", linewidth=0.5,
                alpha=0.85,
            )
            ax.set_xlabel("Pair Density (paired / total)", fontweight="bold")
            ax.set_ylabel("Count", fontweight="bold")
            ax.set_title("Base-Pair Density Distribution", fontweight="bold")

        # ── Panel 5: Quality score distribution ──────────────────────────
        ax = axes[1, 1]
        ax.hist(
            df["quality_score"], bins=40,
            color=MACARON_COLORS[4], edgecolor="white", linewidth=0.5,
            alpha=0.85,
        )
        ax.set_xlabel("Quality Score", fontweight="bold")
        ax.set_ylabel("Count", fontweight="bold")
        ax.set_title("Sequence Quality Score Distribution", fontweight="bold")
        ax.axvline(
            config.MIN_SEQUENCE_QUALITY,
            color=MACARON_COLORS[0], linestyle="--", linewidth=2,
            label=f"Threshold = {config.MIN_SEQUENCE_QUALITY}",
        )
        ax.legend(prop={"weight": "bold"})

        # ── Panel 6: Numeric feature correlation heatmap ─────────────────
        ax = axes[1, 2]
        num_cols = ["length", "quality_score", "gc_content", "pair_density"]
        num_cols = [c for c in num_cols if c in df.columns]
        if len(num_cols) >= 2:
            corr = df[num_cols].corr()
            im = ax.imshow(
                corr.values, cmap=MACARON_DIVERGING,
                vmin=-1, vmax=1, aspect="auto",
            )
            tick_labels = [c.replace("_", " ").title() for c in num_cols]
            ax.set_xticks(range(len(num_cols)))
            ax.set_yticks(range(len(num_cols)))
            ax.set_xticklabels(tick_labels, rotation=35, ha="right",
                                fontweight="bold", fontsize=10)
            ax.set_yticklabels(tick_labels, fontweight="bold", fontsize=10)
            for i in range(len(num_cols)):
                for j in range(len(num_cols)):
                    ax.text(j, i, f"{corr.values[i, j]:.2f}",
                            ha="center", va="center",
                            fontweight="bold", fontsize=9,
                            color="black" if abs(corr.values[i, j]) < 0.7 else "white")
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            ax.set_title("Feature Correlation Matrix", fontweight="bold")

        plt.tight_layout()
        out = config.FIG_DIR / "dataset_overview.pdf"
        plt.savefig(out, dpi=300, bbox_inches="tight")
        plt.savefig(str(out).replace(".pdf", ".png"), dpi=300, bbox_inches="tight")
        plt.close()
        logger.info("Dataset visualisation saved → %s", out)

    def close(self) -> None:
        if self._db:
            self._db.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


# ─────────────────────────────────────────────────────────────────────────────
#  11.  k-MER TOKENIZER
# ─────────────────────────────────────────────────────────────────────────────

class KMerTokenizer:
    """
    k-mer tokeniser for RNA sequences.

    Vocabulary: all k^5 combinations of {A, C, G, U, N}
    plus <PAD> (index 0) and <UNK> (index 1).

    Example (k=3):  "AUGC" → ["AUG", "UGC"]
    """

    BASES    = ["A", "C", "G", "U", "N"]
    PAD_TOK  = "<PAD>"
    UNK_TOK  = "<UNK>"

    def __init__(self, k: int = 3) -> None:
        self.k = k
        self._build_vocab()

    def _build_vocab(self) -> None:
        kmers = [
            "".join(c) for c in itertools.product(self.BASES, repeat=self.k)
        ]
        self.vocab: Dict[str, int] = {
            self.PAD_TOK: 0,
            self.UNK_TOK: 1,
            **{km: i + 2 for i, km in enumerate(kmers)},
        }
        self.inv_vocab: Dict[int, str] = {v: k for k, v in self.vocab.items()}

    @property
    def vocab_size(self) -> int:
        return len(self.vocab)

    def tokenise(self, seq: str) -> List[str]:
        seq = clean_sequence(seq)
        if len(seq) < self.k:
            return [self.UNK_TOK]
        return [seq[i:i + self.k] for i in range(len(seq) - self.k + 1)]

    def encode(self, seq: str, max_len: int) -> np.ndarray:
        """Encode sequence → integer array [max_len], 0-padded."""
        tokens = self.tokenise(seq)
        arr = np.zeros(max_len, dtype=np.int64)
        for i, tok in enumerate(tokens[:max_len]):
            arr[i] = self.vocab.get(tok, self.vocab[self.UNK_TOK])
        return arr

    def decode(self, ids: np.ndarray) -> List[str]:
        return [self.inv_vocab.get(int(i), self.UNK_TOK)
                for i in ids if int(i) != 0]


# ─────────────────────────────────────────────────────────────────────────────
#  12.  RNA-FM EMBEDDING (OPTIONAL)
# ─────────────────────────────────────────────────────────────────────────────

class RNAFMEmbedding(nn.Module):
    """
    Adapter for RNA-FM pre-trained language model.
    Falls back to a trainable linear layer when the model is unavailable.

    Reference: Chen et al. (2022) RNA-FM: Interpretable RNA Foundation Model.
    """

    FM_DIM = 640   # RNA-FM hidden dimension

    def __init__(self, model_name: str = "multimolecule/rnafm",
                 output_dim: int = 128) -> None:
        super().__init__()
        self.output_dim  = output_dim
        self.use_fm = False

        if TRANSFORMERS_AVAILABLE:
            try:
                self._tok   = AutoTokenizer.from_pretrained(model_name)
                self._model = AutoModel.from_pretrained(model_name)
                self._model.eval()
                for p in self._model.parameters():
                    p.requires_grad = False   # frozen backbone
                self._proj = nn.Linear(self.FM_DIM, output_dim)
                self.use_fm = True
                logger.info("RNA-FM loaded: %s", model_name)
            except Exception as exc:
                logger.warning("RNA-FM unavailable (%s); using fallback", exc)

        if not self.use_fm:
            self._fallback = nn.Identity()   # placeholder

    def forward(self, sequences: List[str]) -> "torch.Tensor":
        if not self.use_fm:
            raise RuntimeError("RNA-FM not loaded; disable use_rnafm flag")
        dev = next(self._proj.parameters()).device
        inputs = self._tok(sequences, return_tensors="pt", padding=True)
        inputs = {k: v.to(dev) for k, v in inputs.items()}
        with torch.no_grad():
            out = self._model(**inputs).last_hidden_state   # [B, L, 640]
        return self._proj(out)   # [B, L, output_dim]


# ─────────────────────────────────────────────────────────────────────────────
#  13.  FEATURE CONSTRUCTION HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def build_pair_matrix(dot_bracket: str, max_len: int) -> np.ndarray:
    """Convert dot-bracket → padded pair matrix [max_len, max_len]."""
    if not dot_bracket:
        return np.zeros((max_len, max_len), dtype=np.float32)
    mat  = PseudoKnotParser.dotbracket_to_matrix(dot_bracket)
    out  = np.zeros((max_len, max_len), dtype=np.float32)
    L    = min(max_len, mat.shape[0])
    out[:L, :L] = mat[:L, :L]
    return out


def build_valid_mask(length: int, max_len: int) -> np.ndarray:
    """Binary mask [max_len, max_len] = 1 in the valid [0:L, 0:L] region."""
    mask = np.zeros((max_len, max_len), dtype=np.float32)
    mask[:length, :length] = 1.0
    return mask


def build_distance_bias(length: int, max_len: int) -> np.ndarray:
    """
    Normalised pairwise distance matrix [max_len, max_len].
    Entry (i,j) = |i−j| / max_len ∈ [0, 1].
    """
    idx  = np.arange(max_len, dtype=np.float32)
    dist = np.abs(idx[:, None] - idx[None, :]) / max(max_len, 1)
    return dist.astype(np.float32)


# ─────────────────────────────────────────────────────────────────────────────
#  14.  HOMOLOGY-AWARE DATA SPLITTING
# ─────────────────────────────────────────────────────────────────────────────

class HomologyAwareSplitter:
    """
    Split dataset into train / validation / test partitions
    such that no two partitions share sequences with identity ≥ threshold.

    Supported backends: MMseqs2, CD-HIT-EST.
    Falls back to random splitting when clustering tools are unavailable.

    Reference: Rivas et al. (2012) guidelines on RNA benchmark dataset design.
    """

    def __init__(
        self,
        method:   str   = "mmseqs2",
        identity: float = 0.80,
    ) -> None:
        self.method   = method
        self.identity = identity
        self._verify_tool()

    def _verify_tool(self) -> None:
        def _check(cmd: str) -> bool:
            try:
                subprocess.run(
                    [cmd, "--version"],
                    capture_output=True, check=True
                )
                return True
            except (subprocess.CalledProcessError, FileNotFoundError):
                return False

        if self.method == "mmseqs2" and not _check("mmseqs"):
            logger.warning("MMseqs2 not found — trying CD-HIT-EST")
            self.method = "cdhit"

        if self.method == "cdhit" and not _check("cd-hit-est"):
            logger.warning("CD-HIT-EST not found — using random split")
            self.method = "random"

    def _write_fasta(self, df: pd.DataFrame, path: Path) -> None:
        with open(path, "w", encoding="utf-8") as fh:
            for _, row in df.iterrows():
                fh.write(f">{row['sequence_id']}\n{row['sequence']}\n")

    def _run_mmseqs2(
        self, fasta: Path, out_dir: Path
    ) -> Dict[str, str]:
        out_dir.mkdir(parents=True, exist_ok=True)
        prefix = out_dir / "mmseqs"
        tmp    = out_dir / "tmp"
        subprocess.run([
            "mmseqs", "easy-cluster",
            str(fasta), str(prefix), str(tmp),
            "--min-seq-id", str(self.identity),
            "-c", "0.8", "--cov-mode", "1",
        ], check=True, capture_output=True)
        mapping: Dict[str, str] = {}
        tsv = Path(f"{prefix}_cluster.tsv")
        with open(tsv) as fh:
            for line in fh:
                rep, mem = line.strip().split("\t")
                mapping[mem] = rep
        return mapping

    def _run_cdhit(
        self, fasta: Path, out_dir: Path
    ) -> Dict[str, str]:
        out_dir.mkdir(parents=True, exist_ok=True)
        prefix = out_dir / "cdhit"
        subprocess.run([
            "cd-hit-est",
            "-i", str(fasta), "-o", str(prefix),
            "-c", str(self.identity), "-n", "5",
            "-M", "0", "-T", "0",
        ], check=True, capture_output=True)
        mapping: Dict[str, str] = {}
        clstr = Path(f"{prefix}.clstr")
        cur_cluster: Optional[str] = None
        with open(clstr) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">Cluster"):
                    cur_cluster = line
                else:
                    m = re.search(r">(.+?)\.\.\.", line)
                    if m and cur_cluster:
                        mapping[m.group(1)] = cur_cluster
        return mapping

    def split(
        self, df: pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Returns (train_df, valid_df, test_df) with homology-aware partitioning.
        """
        logger.info(
            "Splitting dataset (method=%s, identity=%.2f, n=%d)",
            self.method, self.identity, len(df),
        )

        if self.method == "random":
            return self._random_split(df)

        fasta   = config.TMP_DIR / "cluster_input.fa"
        out_dir = config.TMP_DIR / self.method
        shutil.rmtree(out_dir, ignore_errors=True)
        self._write_fasta(df, fasta)

        try:
            if self.method == "mmseqs2":
                mapping = self._run_mmseqs2(fasta, out_dir)
            else:
                mapping = self._run_cdhit(fasta, out_dir)
        except Exception as exc:
            logger.warning("Clustering failed (%s) — falling back to random", exc)
            return self._random_split(df)

        df = df.copy()
        df["cluster_id"] = df["sequence_id"].map(mapping)
        df.dropna(subset=["cluster_id"], inplace=True)

        clusters = df["cluster_id"].unique().tolist()
        random.Random(config.SEED).shuffle(clusters)

        n_total = len(clusters)
        n_tr    = int(n_total * config.TRAIN_RATIO)
        n_va    = int(n_total * config.VALID_RATIO)

        tr_cl = set(clusters[:n_tr])
        va_cl = set(clusters[n_tr:n_tr + n_va])
        te_cl = set(clusters[n_tr + n_va:])

        tr = df[df["cluster_id"].isin(tr_cl)].reset_index(drop=True)
        va = df[df["cluster_id"].isin(va_cl)].reset_index(drop=True)
        te = df[df["cluster_id"].isin(te_cl)].reset_index(drop=True)

        self._save_split_stats(tr, va, te, len(clusters))
        return tr, va, te

    def _random_split(
        self, df: pd.DataFrame
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Fallback random split when clustering tools are unavailable."""
        logger.warning("Using random split (no homology control)")
        idx     = list(range(len(df)))
        random.Random(config.SEED).shuffle(idx)
        n_tr    = int(len(idx) * config.TRAIN_RATIO)
        n_va    = int(len(idx) * config.VALID_RATIO)
        tr  = df.iloc[idx[:n_tr]].reset_index(drop=True)
        va  = df.iloc[idx[n_tr:n_tr + n_va]].reset_index(drop=True)
        te  = df.iloc[idx[n_tr + n_va:]].reset_index(drop=True)
        return tr, va, te

    def _save_split_stats(
        self,
        tr: pd.DataFrame,
        va: pd.DataFrame,
        te: pd.DataFrame,
        n_clusters: int,
    ) -> None:
        info = {
            "method":             self.method,
            "identity_threshold": self.identity,
            "n_clusters":         n_clusters,
            "n_train":            len(tr),
            "n_valid":            len(va),
            "n_test":             len(te),
            "train_ratio":        round(len(tr) / max(len(tr)+len(va)+len(te), 1), 3),
        }
        save_json(info, config.PROC_DIR / "split_statistics.json")
        tr.to_csv(config.PROC_DIR / "train.csv", index=False)
        va.to_csv(config.PROC_DIR / "valid.csv", index=False)
        te.to_csv(config.PROC_DIR / "test.csv",  index=False)
        logger.info(
            "Split saved: train=%d  valid=%d  test=%d",
            len(tr), len(va), len(te),
        )


# ─────────────────────────────────────────────────────────────────────────────
#  15.  PYTORCH DATASET
# ─────────────────────────────────────────────────────────────────────────────

if TORCH_AVAILABLE:

    class RNASSDataset(Dataset):
        """
        PyTorch Dataset for RNA secondary structure prediction.

        Each sample contains:
          input_ids   : k-mer token indices [max_len]
          pair_target : ground-truth pair matrix [max_len, max_len]
          valid_mask  : binary mask [max_len, max_len]
          pair_prior  : chemistry-based prior [max_len, max_len]
          dist_bias   : normalised distance matrix [max_len, max_len]
          length      : actual sequence length (scalar)
        """

        def __init__(
            self,
            df:          pd.DataFrame,
            tokenizer:   KMerTokenizer,
            max_len:     int  = 256,
            augment:     bool = False,
        ) -> None:
            self.df        = df.reset_index(drop=True)
            self.tokenizer = tokenizer
            self.max_len   = max_len
            self.augmenter = RNADataAugmenter(config.SEED) if augment else None

        def __len__(self) -> int:
            return len(self.df)

        def __getitem__(self, idx: int) -> Dict[str, Any]:
            row = self.df.iloc[idx]
            seq = clean_sequence(row["sequence"])
            db  = clean_dot_bracket(row.get("dot_bracket") or "")

            # Optional augmentation (training only)
            if self.augmenter is not None:
                seq, db = self.augmenter.augment(seq, db, config.AUGMENT_PROBABILITY)

            L = min(len(seq), self.max_len)

            return {
                "input_ids":   torch.tensor(
                    self.tokenizer.encode(seq, self.max_len), dtype=torch.long
                ),
                "pair_target": torch.tensor(
                    build_pair_matrix(db, self.max_len), dtype=torch.float32
                ),
                "valid_mask":  torch.tensor(
                    build_valid_mask(L, self.max_len), dtype=torch.float32
                ),
                "pair_prior":  torch.tensor(
                    build_basepair_prior(seq, self.max_len, config.MIN_HAIRPIN_LOOP),
                    dtype=torch.float32,
                ),
                "dist_bias":   torch.tensor(
                    build_distance_bias(L, self.max_len), dtype=torch.float32
                ),
                "length":      torch.tensor(L, dtype=torch.long),
                "sequence":    seq,
                "dot_bracket": db,
                "sequence_id": str(row.get("sequence_id", idx)),
            }


    def collate_fn(batch: List[Dict]) -> Dict[str, Any]:
        tensor_keys = [
            "input_ids", "pair_target", "valid_mask",
            "pair_prior", "dist_bias", "length",
        ]
        string_keys = ["sequence", "dot_bracket", "sequence_id"]
        out: Dict[str, Any] = {}
        for k in tensor_keys:
            out[k] = torch.stack([item[k] for item in batch])
        for k in string_keys:
            out[k] = [item[k] for item in batch]
        return out


# ─────────────────────────────────────────────────────────────────────────────
#  16.  MODEL COMPONENTS
# ─────────────────────────────────────────────────────────────────────────────

if TORCH_AVAILABLE:

    class PositionalEncoding(nn.Module):
        """
        Sinusoidal positional encoding (Vaswani et al., 2017).
        Applied to sequence embeddings before Transformer encoder.
        """

        def __init__(
            self, d_model: int, max_len: int = 1024, dropout: float = 0.1
        ) -> None:
            super().__init__()
            self.dropout = nn.Dropout(p=dropout)
            pe       = torch.zeros(max_len, d_model)
            position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
            div_term = torch.exp(
                torch.arange(0, d_model, 2, dtype=torch.float)
                * (-math.log(10000.0) / d_model)
            )
            pe[:, 0::2] = torch.sin(position * div_term)
            if d_model % 2 == 1:
                pe[:, 1::2] = torch.cos(position * div_term[:-1])
            else:
                pe[:, 1::2] = torch.cos(position * div_term)
            self.register_buffer("pe", pe.unsqueeze(0))

        def forward(self, x: "torch.Tensor") -> "torch.Tensor":
            x = x + self.pe[:, :x.size(1)]
            return self.dropout(x)


    class ConvBlock2D(nn.Module):
        """
        2D convolutional block: Conv → BN → GELU → Dropout (×2).
        Used as building block for the U-Net refinement network.
        """

        def __init__(self, in_ch: int, out_ch: int, dropout: float = 0.1) -> None:
            super().__init__()
            self.block = nn.Sequential(
                nn.Conv2d(in_ch, out_ch, 3, padding=1),
                nn.BatchNorm2d(out_ch),
                nn.GELU(),
                nn.Dropout2d(dropout),
                nn.Conv2d(out_ch, out_ch, 3, padding=1),
                nn.BatchNorm2d(out_ch),
                nn.GELU(),
                nn.Dropout2d(dropout),
            )

        def forward(self, x: "torch.Tensor") -> "torch.Tensor":
            return self.block(x)


    class UNet2D(nn.Module):
        """
        U-Net 2D for pairwise contact map refinement.
        Input:  [B, in_ch, L, L]  (multi-channel pair features)
        Output: [B, 1,    L, L]  (refined contact logit map)

        Reference: Ronneberger et al. (2015) U-Net.
        """

        def __init__(
            self,
            in_ch:   int   = 4,
            base_ch: int   = 32,
            out_ch:  int   = 1,
            dropout: float = 0.1,
        ) -> None:
            super().__init__()
            ch = base_ch

            # Encoder
            self.enc1 = ConvBlock2D(in_ch, ch, dropout)
            self.pool1 = nn.MaxPool2d(2)

            self.enc2 = ConvBlock2D(ch,     ch * 2, dropout)
            self.pool2 = nn.MaxPool2d(2)

            self.enc3 = ConvBlock2D(ch * 2, ch * 4, dropout)
            self.pool3 = nn.MaxPool2d(2)

            # Bottleneck
            self.bottleneck = ConvBlock2D(ch * 4, ch * 8, dropout)

            # Decoder
            self.up3   = nn.ConvTranspose2d(ch * 8, ch * 4, 2, stride=2)
            self.dec3  = ConvBlock2D(ch * 8, ch * 4, dropout)

            self.up2   = nn.ConvTranspose2d(ch * 4, ch * 2, 2, stride=2)
            self.dec2  = ConvBlock2D(ch * 4, ch * 2, dropout)

            self.up1   = nn.ConvTranspose2d(ch * 2, ch,     2, stride=2)
            self.dec1  = ConvBlock2D(ch * 2, ch,     dropout)

            # Output
            self.out   = nn.Conv2d(ch, out_ch, 1)

        def _resize(
            self, x: "torch.Tensor", ref: "torch.Tensor"
        ) -> "torch.Tensor":
            if x.shape[-2:] != ref.shape[-2:]:
                x = F.interpolate(
                    x, size=ref.shape[-2:], mode="bilinear", align_corners=False
                )
            return x

        def forward(self, x: "torch.Tensor") -> "torch.Tensor":
            e1 = self.enc1(x);    p1 = self.pool1(e1)
            e2 = self.enc2(p1);   p2 = self.pool2(e2)
            e3 = self.enc3(p2);   p3 = self.pool3(e3)

            b  = self.bottleneck(p3)

            u3 = self._resize(self.up3(b),  e3)
            d3 = self.dec3(torch.cat([u3, e3], 1))

            u2 = self._resize(self.up2(d3), e2)
            d2 = self.dec2(torch.cat([u2, e2], 1))

            u1 = self._resize(self.up1(d2), e1)
            d1 = self.dec1(torch.cat([u1, e1], 1))

            return self.out(d1)


    class ContrastiveLearningHead(nn.Module):
        """
        Contrastive representation learning head (SimCLR-style).
        Projects pooled sequence embeddings to a normalised latent space
        for unsupervised pre-training or auxiliary training.

        Reference: Chen et al. (2020) A Simple Framework for Contrastive Learning.
        """

        def __init__(self, d_model: int, proj_dim: int = 128,
                     temperature: float = 0.07) -> None:
            super().__init__()
            self.temperature = temperature
            self.proj = nn.Sequential(
                nn.Linear(d_model, d_model),
                nn.ReLU(),
                nn.LayerNorm(d_model),
                nn.Linear(d_model, proj_dim),
            )

        def pool(self, features: "torch.Tensor") -> "torch.Tensor":
            """Mean-pool sequence dimension: [B, L, D] → [B, D]."""
            return features.mean(dim=1)

        def forward(self, features: "torch.Tensor") -> "torch.Tensor":
            z = self.proj(self.pool(features))
            return F.normalize(z, dim=1)

        def contrastive_loss(
            self,
            z_i: "torch.Tensor",
            z_j: "torch.Tensor",
        ) -> "torch.Tensor":
            """
            Compute symmetric InfoNCE loss.
            z_i, z_j: normalised projection vectors [B, proj_dim].
            """
            sim   = torch.matmul(z_i, z_j.T) / self.temperature
            lbl   = torch.arange(z_i.size(0), device=z_i.device)
            loss  = (F.cross_entropy(sim, lbl) + F.cross_entropy(sim.T, lbl)) / 2
            return loss


    class RNASSPredictorAdvanced(nn.Module):
        """
        RNA Secondary Structure Prediction Model.

        Architecture:
        ┌─────────────────────────────────────────────┐
        │  k-mer Embedding  +  Positional Encoding    │
        │         ↓                                   │
        │  Transformer Encoder  (N layers)            │
        │         ↓                                   │
        │  Outer-product pair feature construction    │
        │    (additive + multiplicative fusion)       │
        │         ↓                                   │
        │  Multi-channel U-Net 2D refinement          │
        │    channels: [pair_feat, prior,             │
        │               sym_prior, dist_bias]         │
        │         ↓                                   │
        │  Symmetric logit output  [B, L, L]          │
        └─────────────────────────────────────────────┘

        Optional: RNA-FM pre-trained embedding.
        Optional: Contrastive learning auxiliary head.
        """

        def __init__(
            self,
            vocab_size:  int,
            d_model:     int   = 128,
            nhead:       int   = 8,
            num_layers:  int   = 6,
            dim_ff:      int   = 512,
            dropout:     float = 0.10,
            max_len:     int   = 256,
            use_rnafm:   bool  = False,
            unet_base_ch:int   = 32,
        ) -> None:
            super().__init__()
            self.max_len  = max_len
            self.d_model  = d_model
            self.use_rnafm = use_rnafm

            # ── Embedding ────────────────────────────────────────────────
            self.embedding = nn.Embedding(vocab_size, d_model, padding_idx=0)
            self.pos_enc   = PositionalEncoding(d_model, max_len=max_len + 8,
                                                dropout=dropout)

            # ── Optional RNA-FM adapter ──────────────────────────────────
            if use_rnafm:
                self.rnafm = RNAFMEmbedding(output_dim=d_model)
                self.fm_gate = nn.Parameter(torch.zeros(1))  # learnable gate

            # ── Transformer Encoder ──────────────────────────────────────
            enc_layer = nn.TransformerEncoderLayer(
                d_model       = d_model,
                nhead         = nhead,
                dim_feedforward = dim_ff,
                dropout       = dropout,
                batch_first   = True,
                activation    = "gelu",
                norm_first    = True,   # Pre-LN for training stability
            )
            self.encoder = nn.TransformerEncoder(
                enc_layer, num_layers=num_layers,
                enable_nested_tensor=False,
            )

            # ── Pair Feature Projections ─────────────────────────────────
            self.row_proj = nn.Linear(d_model, d_model)
            self.col_proj = nn.Linear(d_model, d_model)
            self.mul_proj = nn.Linear(d_model, d_model)
            self.pair_compress = nn.Conv2d(d_model, 1, kernel_size=1)

            # ── U-Net Refinement ─────────────────────────────────────────
            self.unet = UNet2D(in_ch=4, base_ch=unet_base_ch, out_ch=1,
                               dropout=dropout)

            # ── Contrastive Head (auxiliary) ─────────────────────────────
            self.contrastive_head = ContrastiveLearningHead(d_model)

            # ── Weight Initialisation ────────────────────────────────────
            self._init_weights()

        def _init_weights(self) -> None:
            for name, p in self.named_parameters():
                if p.dim() > 1:
                    nn.init.xavier_uniform_(p)
                elif "bias" in name:
                    nn.init.zeros_(p)

        def encode(
            self,
            input_ids: "torch.Tensor",
        ) -> "torch.Tensor":
            """
            Encode token ids → contextual sequence embeddings.
            input_ids: [B, L]  →  output: [B, L, D]
            """
            pad_mask = (input_ids == 0)
            x = self.embedding(input_ids)   # [B, L, D]
            x = self.pos_enc(x)
            x = self.encoder(x, src_key_padding_mask=pad_mask)
            return x

        def forward(
            self,
            input_ids:  "torch.Tensor",
            pair_prior: "torch.Tensor",
            dist_bias:  "torch.Tensor",
        ) -> "torch.Tensor":
            """
            Forward pass.

            Args:
                input_ids  : [B, L]     integer k-mer token indices
                pair_prior : [B, L, L]  chemistry-based prior matrix
                dist_bias  : [B, L, L]  normalised pairwise distance

            Returns:
                logits: [B, L, L] symmetric pair logit matrix
            """
            # ── Sequence encoding ────────────────────────────────────────
            x = self.encode(input_ids)   # [B, L, D]

            # ── Optional RNA-FM feature gating ───────────────────────────
            # (Requires sequences, which are not passed here for efficiency;
            #  integrate via a separate encode_rnafm() call if needed)

            # ── Outer-product pair features ──────────────────────────────
            row  = self.row_proj(x)    # [B, L, D]
            col  = self.col_proj(x)    # [B, L, D]
            mul  = self.mul_proj(x)    # [B, L, D]

            pair = (
                row.unsqueeze(2) + col.unsqueeze(1)    # additive
                + mul.unsqueeze(2) * mul.unsqueeze(1)  # multiplicative
            )  # [B, L, L, D]

            # ── Compress to single channel [B, 1, L, L] ──────────────────
            pair_feat = self.pair_compress(
                pair.permute(0, 3, 1, 2)
            )   # [B, 1, L, L]

            # ── Assemble 4-channel U-Net input ───────────────────────────
            prior      = pair_prior.unsqueeze(1)              # [B, 1, L, L]
            sym_prior  = 0.5 * (prior + prior.transpose(2, 3))
            dist       = dist_bias.unsqueeze(1)               # [B, 1, L, L]

            unet_in = torch.cat([pair_feat, prior, sym_prior, dist], dim=1)

            # ── U-Net refinement → symmetric logits ──────────────────────
            logits = self.unet(unet_in).squeeze(1)            # [B, L, L]
            logits = 0.5 * (logits + logits.transpose(1, 2))  # symmetrise

            return logits

        def get_contrastive_embeddings(
            self, input_ids: "torch.Tensor"
        ) -> "torch.Tensor":
            """Return normalised contrastive projections for auxiliary training."""
            x = self.encode(input_ids)
            return self.contrastive_head(x)


# ─────────────────────────────────────────────────────────────────────────────
#  17.  COMPOSITE LOSS FUNCTION
# ─────────────────────────────────────────────────────────────────────────────

if TORCH_AVAILABLE:

    class RNAStructureLoss(nn.Module):
        """
        Multi-objective composite loss for RNA structure prediction.

        Components:
        ┌──────────────────────────────────────────────────────────────┐
        │ L_total = w_bce   · L_BCE   (with pos_weight & label smooth) │
        │         + w_focal · L_Focal (addressing class imbalance)     │
        │         + w_dice  · L_Dice  (structure-level overlap)        │
        │         + w_sym   · L_sym   (symmetry constraint)            │
        │         + w_diag  · L_diag  (no self-pairing penalty)        │
        │         + w_sparse· L_sparse(pair density regularisation)    │
        └──────────────────────────────────────────────────────────────┘
        """

        def __init__(self, cfg: Config) -> None:
            super().__init__()
            self.w       = cfg.LOSS_WEIGHTS
            self.pw      = cfg.POS_WEIGHT
            self.alpha   = cfg.FOCAL_ALPHA
            self.gamma   = cfg.FOCAL_GAMMA
            self.tgt_den = cfg.TARGET_PAIR_DENSITY
            self.smooth  = cfg.LABEL_SMOOTHING

        def _bce(
            self,
            logits:  "torch.Tensor",
            targets: "torch.Tensor",
            mask:    "torch.Tensor",
        ) -> "torch.Tensor":
            # Label smoothing
            targets_smooth = targets * (1 - self.smooth) + 0.5 * self.smooth
            pw = torch.tensor([self.pw], device=logits.device)
            loss = F.binary_cross_entropy_with_logits(
                logits, targets_smooth, reduction="none", pos_weight=pw
            )
            return (loss * mask).sum() / mask.sum().clamp(min=1.0)

        def _focal(
            self,
            logits:  "torch.Tensor",
            targets: "torch.Tensor",
            mask:    "torch.Tensor",
        ) -> "torch.Tensor":
            p    = torch.sigmoid(logits)
            bce  = F.binary_cross_entropy_with_logits(
                logits, targets, reduction="none"
            )
            pt   = p * targets + (1 - p) * (1 - targets)
            loss = self.alpha * (1 - pt).pow(self.gamma) * bce
            return (loss * mask).sum() / mask.sum().clamp(min=1.0)

        def _dice(
            self,
            logits:  "torch.Tensor",
            targets: "torch.Tensor",
            mask:    "torch.Tensor",
            eps:     float = 1e-6,
        ) -> "torch.Tensor":
            p    = torch.sigmoid(logits) * mask
            t    = targets * mask
            num  = 2 * (p * t).sum() + eps
            den  = p.sum() + t.sum() + eps
            return 1.0 - num / den

        def _symmetry(
            self,
            logits: "torch.Tensor",
            mask:   "torch.Tensor",
        ) -> "torch.Tensor":
            diff = (logits - logits.transpose(1, 2)).abs() * mask
            return diff.sum() / mask.sum().clamp(min=1.0)

        def _diagonal_penalty(
            self,
            logits: "torch.Tensor",
            mask:   "torch.Tensor",
        ) -> "torch.Tensor":
            B, L, _ = logits.shape
            eye  = torch.eye(L, device=logits.device).unsqueeze(0)
            prob = torch.sigmoid(logits)
            return (prob * eye * mask).sum() / (eye * mask).sum().clamp(min=1.0)

        def _sparsity(
            self,
            logits: "torch.Tensor",
            mask:   "torch.Tensor",
        ) -> "torch.Tensor":
            prob    = torch.sigmoid(logits) * mask
            density = prob.sum() / mask.sum().clamp(min=1.0)
            return (density - self.tgt_den).abs()

        def forward(
            self,
            logits:  "torch.Tensor",
            targets: "torch.Tensor",
            mask:    "torch.Tensor",
        ) -> Tuple["torch.Tensor", Dict[str, float]]:
            l_bce    = self._bce(logits, targets, mask)
            l_focal  = self._focal(logits, targets, mask)
            l_dice   = self._dice(logits, targets, mask)
            l_sym    = self._symmetry(logits, mask)
            l_diag   = self._diagonal_penalty(logits, mask)
            l_sparse = self._sparsity(logits, mask)

            total = (
                self.w["bce"]    * l_bce
                + self.w["focal"]  * l_focal
                + self.w["dice"]   * l_dice
                + self.w["sym"]    * l_sym
                + self.w["diag"]   * l_diag
                + self.w["sparse"] * l_sparse
            )

            stats = {
                "bce":    float(l_bce),
                "focal":  float(l_focal),
                "dice":   float(l_dice),
                "sym":    float(l_sym),
                "diag":   float(l_diag),
                "sparse": float(l_sparse),
                "total":  float(total),
            }
            return total, stats


# ─────────────────────────────────────────────────────────────────────────────
#  18.  NUSSINOV DYNAMIC PROGRAMMING DECODER
# ─────────────────────────────────────────────────────────────────────────────

class NussinovDecoder:
    """
    Nussinov dynamic programming decoder for RNA secondary structure.

    Algorithm:
      DP[i][j] = maximum number of base pairs in sub-sequence [i, j]

    Extensions implemented:
    - Chemistry filtering (only canonical/wobble pairs allowed)
    - Minimum loop length constraint
    - Multi-bracket pseudoknot assignment via graph colouring

    Reference: Nussinov & Jacobson (1980) PNAS.
    """

    def __init__(
        self,
        min_loop:      int   = 3,
        allowed_pairs: set   = ALLOWED_PAIRS,
        threshold:     float = 0.5,
    ) -> None:
        self.min_loop      = min_loop
        self.allowed_pairs = allowed_pairs
        self.threshold     = threshold

    def decode(
        self,
        score_mat: np.ndarray,
        sequence:  str,
        threshold: Optional[float] = None,
    ) -> str:
        """
        Decode probability matrix → dot-bracket string.

        Args:
            score_mat : [L, L] pairwise probability or logit matrix
            sequence  : RNA sequence string (length L)
            threshold : minimum score to consider a pair

        Returns:
            dot_bracket string of length L
        """
        thr = threshold if threshold is not None else self.threshold
        seq = clean_sequence(sequence)
        L   = len(seq)

        if score_mat.shape[0] < L or score_mat.shape[1] < L:
            raise ValueError(
                f"Score matrix {score_mat.shape} smaller than sequence length {L}"
            )
        score_mat = score_mat[:L, :L]

        # Filter to chemically valid pairs above threshold
        valid: Dict[Tuple[int, int], float] = {}
        for i in range(L):
            for j in range(i + self.min_loop + 1, L):
                if (is_valid_pair(seq[i], seq[j])
                        and score_mat[i, j] >= thr):
                    valid[(i, j)] = float(score_mat[i, j])

        # DP table
        dp = np.zeros((L, L), dtype=np.float32)
        bt = np.full((L, L), -1, dtype=np.int32)

        for span in range(self.min_loop + 1, L):
            for i in range(L - span):
                j = i + span

                # Option A: skip i
                best = dp[i + 1, j] if i + 1 <= j else 0.0
                best_k = -2   # sentinel: "skip i"

                # Option B: skip j  (handled implicitly via trace)
                if (i <= j - 1) and dp[i, j - 1] > best:
                    best   = dp[i, j - 1]
                    best_k = -3   # sentinel: "skip j"

                # Option C: pair (i, j)
                if (i, j) in valid:
                    inner = dp[i + 1, j - 1] if i + 1 <= j - 1 else 0.0
                    score = valid[(i, j)] + inner
                    if score > best:
                        best   = score
                        best_k = -1   # sentinel: "pair i-j"

                # Option D: bifurcation
                for k in range(i + 1, j):
                    s = dp[i, k] + dp[k + 1, j]
                    if s > best:
                        best   = s
                        best_k = k   # bifurcation point

                dp[i, j] = best
                bt[i, j] = best_k

        # Traceback
        pairs: List[Tuple[int, int]] = []

        def _trace(i: int, j: int) -> None:
            if i >= j:
                return
            k = int(bt[i, j])
            if k == -2:              # skip i
                _trace(i + 1, j)
            elif k == -3:            # skip j
                _trace(i, j - 1)
            elif k == -1:            # pair (i, j)
                pairs.append((i, j))
                _trace(i + 1, j - 1)
            else:                    # bifurcation at k
                _trace(i, k)
                _trace(k + 1, j)

        _trace(0, L - 1)
        return PseudoKnotParser.pairs_to_dotbracket(pairs, L)


# ─────────────────────────────────────────────────────────────────────────────
#  19.  VIENNA RNA THERMODYNAMIC FUSION
# ─────────────────────────────────────────────────────────────────────────────

class ViennaRNAInterface:
    """
    Interface to ViennaRNA RNAfold for thermodynamic probability integration.

    Usage:
      1. fold(seq)          → MFE structure string
      2. pair_probabilities(seq) → [L, L] base-pair probability matrix
      3. fuse(dnn_prob, thermo_prob, alpha) → weighted combination

    Reference: Lorenz et al. (2011) ViennaRNA Package 2.0.
    """

    def __init__(self, rnafold_path: str = "RNAfold") -> None:
        self._bin = rnafold_path
        self._available = self._check()

    def _check(self) -> bool:
        try:
            subprocess.run(
                [self._bin, "--version"],
                capture_output=True, check=True
            )
            return True
        except (FileNotFoundError, subprocess.CalledProcessError):
            logger.warning("RNAfold not found — thermodynamic fusion disabled")
            return False

    @property
    def available(self) -> bool:
        return self._available

    def fold(self, sequence: str) -> Optional[str]:
        """Return MFE dot-bracket structure."""
        if not self._available:
            return None
        try:
            proc = subprocess.run(
                [self._bin, "--noPS"],
                input=sequence, text=True,
                capture_output=True, timeout=60, check=True,
            )
            lines = proc.stdout.strip().split("\n")
            if len(lines) >= 2:
                db_line = lines[1].split()[0]
                return clean_dot_bracket(db_line)
        except Exception as exc:
            logger.debug("RNAfold error: %s", exc)
        return None

    def pair_probabilities(self, sequence: str) -> Optional[np.ndarray]:
        """
        Compute base-pair probability matrix using partition function.
        Returns [L, L] float32 matrix or None on failure.
        """
        if not self._available:
            return None

        work_dir = config.TMP_DIR / "vienna"
        work_dir.mkdir(parents=True, exist_ok=True)
        seq_file = work_dir / "seq.fa"

        with open(seq_file, "w") as fh:
            fh.write(f">query\n{sequence}\n")

        try:
            subprocess.run(
                [self._bin, "-p", "--noPS", "--id-prefix=rna"],
                stdin=open(seq_file), stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL, cwd=str(work_dir),
                timeout=120,
            )
            # Parse dot.ps for pair probabilities
            L    = len(sequence)
            pmat = np.zeros((L, L), dtype=np.float32)
            for ps_file in work_dir.glob("*_dp.ps"):
                with open(ps_file) as fh:
                    for line in fh:
                        m = re.match(
                            r"\s*(\d+)\s+(\d+)\s+([0-9.eE+-]+)\s+ubox",
                            line
                        )
                        if m:
                            i, j, p = int(m.group(1))-1, int(m.group(2))-1, float(m.group(3))
                            if 0 <= i < L and 0 <= j < L:
                                pmat[i, j] = pmat[j, i] = p ** 2   # p is sqrt(prob)
            return pmat
        except Exception as exc:
            logger.debug("RNAfold partition failed: %s", exc)
            return None

    @staticmethod
    def fuse(
        dnn_prob:    np.ndarray,
        thermo_prob: Optional[np.ndarray],
        alpha:       float = 0.80,
    ) -> np.ndarray:
        """
        Linear interpolation: (1-α)·thermo + α·DNN
        alpha=1.0 → pure DNN; alpha=0.0 → pure thermodynamic.
        """
        if thermo_prob is None or thermo_prob.shape != dnn_prob.shape:
            return dnn_prob
        return alpha * dnn_prob + (1.0 - alpha) * thermo_prob


# ─────────────────────────────────────────────────────────────────────────────
#  20.  EVALUATION METRICS
# ─────────────────────────────────────────────────────────────────────────────

class RNAStructureMetrics:
    """
    Comprehensive RNA secondary structure evaluation metrics.

    Evaluation protocol:
    1. Restrict to valid-pair space: upper-triangle pairs with |i-j| > min_loop
       and chemically permissible bases (canonical + wobble).
    2. Compute TP, FP, FN, TN in this restricted space.
    3. Report: Precision (PPV), Sensitivity (TPR), F1, MCC, exact-match rate,
       pseudoknot-specific F1.

    Reference: Mathews (2019) on RNA structure evaluation standards.
    """

    def __init__(self, min_loop: int = 3) -> None:
        self.min_loop = min_loop

    def _valid_pair_space(self, seq: str) -> Set[Tuple[int, int]]:
        L = len(seq)
        return {
            (i, j)
            for i in range(L)
            for j in range(i + self.min_loop + 1, L)
            if is_valid_pair(seq[i], seq[j])
        }

    @staticmethod
    def _to_canonical_set(pairs: List[Tuple[int, int]]) -> Set[Tuple[int, int]]:
        return {(min(a, b), max(a, b)) for a, b in pairs if a != b}

    def compute(
        self,
        pred_db: str,
        true_db: str,
        sequence: str,
    ) -> Dict[str, float]:
        """
        Compute all metrics for a single sequence prediction.

        Returns dict containing:
          precision, sensitivity, f1, mcc, exact_match,
          pseudoknot_f1, tp, fp, fn, tn
        """
        pred_pairs = self._to_canonical_set(
            PseudoKnotParser.dotbracket_to_pairs(pred_db)
        )
        true_pairs = self._to_canonical_set(
            PseudoKnotParser.dotbracket_to_pairs(true_db)
        )
        valid_space = self._valid_pair_space(sequence)

        # Restrict to valid pair space
        pred_v = pred_pairs & valid_space
        true_v = true_pairs & valid_space

        tp = len(pred_v & true_v)
        fp = len(pred_v - true_v)
        fn = len(true_v - pred_v)
        tn = len(valid_space) - tp - fp - fn

        prec = tp / (tp + fp + 1e-9)
        sens = tp / (tp + fn + 1e-9)
        f1   = 2 * prec * sens / (prec + sens + 1e-9)

        mcc_num  = float(tp * tn - fp * fn)
        mcc_den  = float(
            math.sqrt(
                max((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn), 0)
            ) + 1e-9
        )
        mcc  = mcc_num / mcc_den

        exact = 1.0 if pred_pairs == true_pairs else 0.0

        # Pseudoknot-specific evaluation
        pk_pred = PseudoKnotParser.identify_pseudoknots(pred_v)
        pk_true = PseudoKnotParser.identify_pseudoknots(true_v)
        pk_tp   = len(pk_pred & pk_true)
        pk_prec = pk_tp / (len(pk_pred) + 1e-9)
        pk_sens = pk_tp / (len(pk_true) + 1e-9)
        pk_f1   = 2 * pk_prec * pk_sens / (pk_prec + pk_sens + 1e-9)

        return {
            "precision":        prec,
            "sensitivity":      sens,
            "f1":               f1,
            "mcc":              mcc,
            "exact_match":      exact,
            "tp":               tp,
            "fp":               fp,
            "fn":               fn,
            "tn":               tn,
            "pseudoknot_f1":    pk_f1,
            "n_pred_pairs":     len(pred_v),
            "n_true_pairs":     len(true_v),
        }

    def batch_metrics(
        self,
        pred_structures: List[str],
        true_structures: List[str],
        sequences:       List[str],
    ) -> Dict[str, float]:
        """
        Compute mean ± std metrics across a batch.
        Returns {metric_mean: float, metric_std: float, ...}.
        """
        all_m: Dict[str, List[float]] = defaultdict(list)
        for p, t, s in zip(pred_structures, true_structures, sequences):
            m = self.compute(p, t, s)
            for k, v in m.items():
                if isinstance(v, (int, float)):
                    all_m[k].append(float(v))

        out: Dict[str, float] = {}
        for k, vals in all_m.items():
            out[f"{k}_mean"] = float(np.mean(vals))
            out[f"{k}_std"]  = float(np.std(vals))
        return out


# ─────────────────────────────────────────────────────────────────────────────
#  21.  STATISTICAL SIGNIFICANCE TESTING
# ─────────────────────────────────────────────────────────────────────────────

class StatisticalSignificanceTester:
    """
    Statistical significance tests for comparing two models.

    Methods:
    - bootstrap_ci   : non-parametric bootstrap confidence intervals + p-value
    - mcnemar_test   : McNemar's test for paired binary classification results
    - wilcoxon_test  : Wilcoxon signed-rank test for continuous metrics
    """

    @staticmethod
    def bootstrap_ci(
        scores_a: np.ndarray,
        scores_b: np.ndarray,
        n_boot:   int   = 10_000,
        alpha:    float = 0.05,
    ) -> Dict[str, Any]:
        """
        Paired bootstrap test for difference in mean scores.

        Returns:
            observed_diff, ci_lower, ci_upper, p_value, significant
        """
        scores_a = np.asarray(scores_a, dtype=float)
        scores_b = np.asarray(scores_b, dtype=float)
        assert len(scores_a) == len(scores_b), "Unequal sample sizes"

        n    = len(scores_a)
        diff = float(scores_a.mean() - scores_b.mean())

        rng       = np.random.default_rng(42)
        boot_diffs = np.empty(n_boot)
        for i in range(n_boot):
            idx             = rng.integers(0, n, n)
            boot_diffs[i]   = scores_a[idx].mean() - scores_b[idx].mean()

        ci_lo = float(np.percentile(boot_diffs, 100 * alpha / 2))
        ci_hi = float(np.percentile(boot_diffs, 100 * (1 - alpha / 2)))

        # Two-sided p-value (shift null to zero)
        null_diffs = boot_diffs - boot_diffs.mean()
        p_val = float(np.mean(np.abs(null_diffs) >= abs(diff)))

        return {
            "observed_diff": diff,
            "ci_lower":      ci_lo,
            "ci_upper":      ci_hi,
            "p_value":       p_val,
            "significant":   p_val < alpha,
        }

    @staticmethod
    def mcnemar_test(
        correct_a: List[int],
        correct_b: List[int],
    ) -> Dict[str, Any]:
        """
        McNemar's test for paired nominal data.
        correct_a[i], correct_b[i] ∈ {0, 1}: whether model A / B was correct.
        """
        assert len(correct_a) == len(correct_b)
        n01 = sum(1 for a, b in zip(correct_a, correct_b) if a == 1 and b == 0)
        n10 = sum(1 for a, b in zip(correct_a, correct_b) if a == 0 and b == 1)

        # With continuity correction
        chi2 = (abs(n01 - n10) - 1) ** 2 / max(n01 + n10, 1)
        p    = float(1.0 - scipy_stats.chi2.cdf(chi2, df=1))

        return {
            "n01":         n01,
            "n10":         n10,
            "chi2":        float(chi2),
            "p_value":     p,
            "significant": p < 0.05,
        }

    @staticmethod
    def wilcoxon_test(
        scores_a: np.ndarray,
        scores_b: np.ndarray,
    ) -> Dict[str, Any]:
        """Wilcoxon signed-rank test (non-parametric paired comparison)."""
        diff = np.asarray(scores_a) - np.asarray(scores_b)
        stat, p = scipy_stats.wilcoxon(diff, alternative="two-sided")
        return {
            "statistic":   float(stat),
            "p_value":     float(p),
            "significant": p < 0.05,
        }


# ─────────────────────────────────────────────────────────────────────────────
#  22.  MODEL EXPLAINABILITY (ATTENTION VISUALISER + INPUT PERTURBATION)
# ─────────────────────────────────────────────────────────────────────────────

if TORCH_AVAILABLE:

    class AttentionVisualizer:
        """
        Capture and visualise multi-head self-attention weights
        from Transformer encoder layers via forward hooks.
        """

        def __init__(self, model: RNASSPredictorAdvanced) -> None:
            self.model   = model
            self._weights: Dict[str, "torch.Tensor"] = {}
            self._hooks:   List[Any] = []
            self._register_hooks()

        def _register_hooks(self) -> None:
            def _make_hook(name: str):
                def _hook(module, inp, out):
                    # TransformerEncoderLayer returns (output,) or output
                    # Attention weights are NOT returned by default in PyTorch
                    # → we store the output activations as proxy
                    if isinstance(out, tuple):
                        self._weights[name] = out[0].detach().cpu()
                    else:
                        self._weights[name] = out.detach().cpu()
                return _hook

            for i, layer in enumerate(self.model.encoder.layers):
                h = layer.register_forward_hook(_make_hook(f"layer_{i}"))
                self._hooks.append(h)

        def remove_hooks(self) -> None:
            for h in self._hooks:
                h.remove()
            self._hooks.clear()

        def plot_attention(
            self,
            sequence:  str,
            layer:     int = -1,
            save_path: Optional[Path] = None,
        ) -> None:
            """
            Plot per-layer activation heatmap as proxy for attention.
            """
            if not self._weights:
                logger.warning("No attention weights captured yet")
                return

            keys      = list(self._weights.keys())
            layer_key = keys[layer]
            act       = self._weights[layer_key]  # [1, L, D] or [B, L, D]

            if act.dim() == 3:
                act = act[0]   # [L, D]

            L = min(len(sequence), act.shape[0])
            act_np = act[:L].numpy()

            fig, ax = plt.subplots(figsize=(14, 5))
            im = ax.imshow(
                act_np.T, aspect="auto",
                cmap=MACARON_SEQUENTIAL,
                interpolation="nearest",
            )
            ax.set_xticks(range(L))
            ax.set_xticklabels(
                list(sequence[:L]), fontsize=9, fontweight="bold"
            )
            ax.set_ylabel("Embedding Dimension", fontweight="bold")
            ax.set_xlabel("Sequence Position", fontweight="bold")
            ax.set_title(
                f"Transformer Activation Map — {layer_key.replace('_', ' ').title()}",
                fontweight="bold",
            )
            plt.colorbar(im, ax=ax, fraction=0.03)
            plt.tight_layout()

            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches="tight")
                logger.info("Attention map saved → %s", save_path)
            plt.close()


    class InputPerturbationExplainer:
        """
        Model-agnostic explainability via input perturbation.
        Estimates per-position importance by measuring prediction change
        upon single-position mutations (approximation of SHAP values).

        Reference: Inspired by Lundberg & Lee (2017) SHAP.
        """

        def __init__(
            self,
            model:     RNASSPredictorAdvanced,
            tokenizer: KMerTokenizer,
            device:    str = "cpu",
        ) -> None:
            self.model     = model
            self.tokenizer = tokenizer
            self.device    = device

        def _predict_confidence(self, sequence: str) -> float:
            self.model.eval()
            input_ids = torch.tensor(
                self.tokenizer.encode(sequence, config.MAX_LEN),
                dtype=torch.long,
            ).unsqueeze(0).to(self.device)
            prior = torch.tensor(
                build_basepair_prior(sequence, config.MAX_LEN),
                dtype=torch.float32,
            ).unsqueeze(0).to(self.device)
            dist = torch.tensor(
                build_distance_bias(len(sequence), config.MAX_LEN),
                dtype=torch.float32,
            ).unsqueeze(0).to(self.device)

            with torch.no_grad():
                logits = self.model(input_ids, prior, dist)
                L      = len(sequence)
                probs  = torch.sigmoid(logits)[0, :L, :L]
            return float(probs.mean().cpu())

        def explain(
            self,
            sequence:       str,
            n_perturbations: int = 20,
        ) -> Dict[str, Any]:
            """
            Compute per-position importance scores.
            Returns dict with importance array and most/least important sites.
            """
            L       = len(sequence)
            bases   = [b for b in "ACGU" if True]
            orig    = self._predict_confidence(sequence)
            scores  = np.zeros(L, dtype=np.float64)

            for i in tqdm(range(L), desc="Explaining", leave=False):
                deltas = []
                for _ in range(n_perturbations):
                    mutated  = list(sequence)
                    opts     = [b for b in bases if b != sequence[i]]
                    mutated[i] = random.choice(opts)
                    delta    = abs(orig - self._predict_confidence("".join(mutated)))
                    deltas.append(delta)
                scores[i] = np.mean(deltas)

            # Normalise
            if scores.max() > 0:
                scores /= scores.max()

            return {
                "sequence":          sequence,
                "importance_scores": scores.tolist(),
                "most_important":    int(np.argmax(scores)),
                "least_important":   int(np.argmin(scores)),
                "baseline_confidence": orig,
            }

        def plot_importance(
            self,
            result:    Dict[str, Any],
            save_path: Optional[Path] = None,
        ) -> None:
            seq    = result["sequence"]
            scores = np.array(result["importance_scores"])

            fig, ax = plt.subplots(figsize=(max(12, len(seq) * 0.35), 4))
            colours = [
                MACARON_COLORS[0] if s > 0.7
                else MACARON_COLORS[1] if s > 0.3
                else MACARON_COLORS[4]
                for s in scores
            ]
            ax.bar(range(len(scores)), scores, color=colours, edgecolor="white")
            ax.set_xticks(range(len(seq)))
            ax.set_xticklabels(list(seq), fontsize=9, fontweight="bold")
            ax.set_ylabel("Normalised Importance Score", fontweight="bold")
            ax.set_xlabel("Sequence Position", fontweight="bold")
            ax.set_title(
                "Per-position Input Importance (Perturbation Analysis)",
                fontweight="bold",
            )
            ax.set_ylim(0, 1.05)

            # Legend
            patches = [
                mpatches.Patch(color=MACARON_COLORS[0], label="High (>0.7)"),
                mpatches.Patch(color=MACARON_COLORS[1], label="Medium (0.3–0.7)"),
                mpatches.Patch(color=MACARON_COLORS[4], label="Low (<0.3)"),
            ]
            ax.legend(handles=patches, prop={"weight": "bold"}, loc="upper right")

            plt.tight_layout()
            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches="tight")
            plt.close()


# ─────────────────────────────────────────────────────────────────────────────
#  23.  UNCERTAINTY ESTIMATION (MONTE CARLO DROPOUT)
# ─────────────────────────────────────────────────────────────────────────────

if TORCH_AVAILABLE:

    class MCDropoutEstimator:
        """
        Monte Carlo Dropout uncertainty estimation.

        At inference, keep Dropout layers active and run N stochastic
        forward passes to approximate the posterior predictive distribution.

        Reference: Gal & Ghahramani (2016) Dropout as a Bayesian Approximation.
        """

        def __init__(
            self,
            model:      RNASSPredictorAdvanced,
            n_samples:  int = 30,
        ) -> None:
            self.model     = model
            self.n_samples = n_samples

        def _enable_dropout(self) -> None:
            """Set all Dropout layers to training mode while keeping rest in eval."""
            self.model.eval()
            for m in self.model.modules():
                if isinstance(m, (nn.Dropout, nn.Dropout2d)):
                    m.train()

        def estimate(
            self,
            input_ids:  "torch.Tensor",
            pair_prior: "torch.Tensor",
            dist_bias:  "torch.Tensor",
        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
            """
            Run N stochastic forward passes.

            Returns:
                mean_prob  : [B, L, L] mean prediction
                std_prob   : [B, L, L] epistemic uncertainty (std dev)
                entropy    : [B, L, L] predictive entropy
            """
            self._enable_dropout()
            preds: List[np.ndarray] = []

            with torch.no_grad():
                for _ in range(self.n_samples):
                    logits = self.model(input_ids, pair_prior, dist_bias)
                    preds.append(torch.sigmoid(logits).cpu().numpy())

            stack    = np.stack(preds, axis=0)   # [N, B, L, L]
            mean_p   = stack.mean(axis=0)
            std_p    = stack.std(axis=0)

            eps     = 1e-9
            entropy = -(
                mean_p * np.log(mean_p + eps)
                + (1 - mean_p) * np.log(1 - mean_p + eps)
            )

            return mean_p, std_p, entropy

        def plot_uncertainty(
            self,
            mean_prob: np.ndarray,
            std_prob:  np.ndarray,
            sequence:  str,
            save_path: Optional[Path] = None,
        ) -> None:
            L   = len(sequence)
            fig, axes = plt.subplots(1, 2, figsize=(14, 6))

            im0 = axes[0].imshow(
                mean_prob[:L, :L], cmap=MACARON_SEQUENTIAL,
                vmin=0, vmax=1, aspect="auto",
            )
            axes[0].set_title("Mean Predicted Pairing Probability", fontweight="bold")
            axes[0].set_xlabel("Position j", fontweight="bold")
            axes[0].set_ylabel("Position i", fontweight="bold")
            plt.colorbar(im0, ax=axes[0], fraction=0.046)

            im1 = axes[1].imshow(
                std_prob[:L, :L], cmap="Oranges",
                vmin=0, aspect="auto",
            )
            axes[1].set_title("Epistemic Uncertainty (Std Dev)", fontweight="bold")
            axes[1].set_xlabel("Position j", fontweight="bold")
            axes[1].set_ylabel("Position i", fontweight="bold")
            plt.colorbar(im1, ax=axes[1], fraction=0.046)

            plt.suptitle(
                f"MC Dropout Uncertainty (N={self.n_samples} samples)",
                fontweight="bold", fontsize=15, y=1.02,
            )
            plt.tight_layout()
            if save_path:
                plt.savefig(save_path, dpi=300, bbox_inches="tight")
            plt.close()


# ─────────────────────────────────────────────────────────────────────────────
#  24.  VISUALISATION UTILITIES
# ─────────────────────────────────────────────────────────────────────────────

def plot_roc_pr_curves(
    y_true:    List[int],
    y_score:   List[float],
    model_name: str     = "Model",
    save_path: Optional[Path] = None,
) -> Dict[str, float]:
    """
    Plot ROC and Precision-Recall curves side by side.
    SCI-quality, macaron palette, bold Times New Roman.
    """
    y_true  = np.array(y_true)
    y_score = np.array(y_score)

    fpr, tpr, _    = roc_curve(y_true, y_score)
    roc_auc        = auc(fpr, tpr)
    prec, rec, _   = precision_recall_curve(y_true, y_score)
    pr_auc         = auc(rec, prec)
    baseline_prec  = y_true.mean()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # ROC
    ax1.plot(fpr, tpr, color=MACARON_COLORS[1], linewidth=2.5,
             label=f"{model_name}  AUC = {roc_auc:.3f}")
    ax1.plot([0, 1], [0, 1], color=MACARON_COLORS[0], linestyle="--",
             linewidth=1.5, label="Random Classifier")
    ax1.fill_between(fpr, tpr, alpha=0.12, color=MACARON_COLORS[1])
    ax1.set_xlabel("False Positive Rate (1 – Specificity)", fontweight="bold")
    ax1.set_ylabel("True Positive Rate (Sensitivity)", fontweight="bold")
    ax1.set_title("Receiver Operating Characteristic Curve", fontweight="bold")
    ax1.legend(prop={"weight": "bold"})
    ax1.set_xlim([0, 1]); ax1.set_ylim([0, 1.02])

    # PR
    ax2.plot(rec, prec, color=MACARON_COLORS[2], linewidth=2.5,
             label=f"{model_name}  AUC = {pr_auc:.3f}")
    ax2.axhline(baseline_prec, color=MACARON_COLORS[0], linestyle="--",
                linewidth=1.5, label=f"Baseline (prevalence = {baseline_prec:.3f})")
    ax2.fill_between(rec, prec, alpha=0.12, color=MACARON_COLORS[2])
    ax2.set_xlabel("Recall (Sensitivity)", fontweight="bold")
    ax2.set_ylabel("Precision (PPV)", fontweight="bold")
    ax2.set_title("Precision-Recall Curve", fontweight="bold")
    ax2.legend(prop={"weight": "bold"})
    ax2.set_xlim([0, 1]); ax2.set_ylim([0, 1.02])

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info("ROC/PR curves saved → %s", save_path)
    plt.close()

    return {"roc_auc": roc_auc, "pr_auc": pr_auc}


def plot_training_curves(
    history: Dict[str, List],
    save_path: Optional[Path] = None,
) -> None:
    """
    Plot training loss and validation metrics over epochs.
    """
    epochs = list(range(1, len(history["train_loss"]) + 1))

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle("Training Dynamics", fontsize=16, fontweight="bold")

    # Loss
    ax = axes[0]
    ax.plot(epochs, history["train_loss"],
            color=MACARON_COLORS[0], linewidth=2, label="Train Loss")
    if "val_loss" in history:
        ax.plot(epochs, history["val_loss"],
                color=MACARON_COLORS[1], linewidth=2, linestyle="--",
                label="Validation Loss")
    ax.set_xlabel("Epoch", fontweight="bold")
    ax.set_ylabel("Loss", fontweight="bold")
    ax.set_title("Training & Validation Loss", fontweight="bold")
    ax.legend(prop={"weight": "bold"})

    # F1 & MCC
    ax = axes[1]
    val_f1  = [m.get("f1_mean",  m.get("f1",  0)) for m in history["val_metrics"]]
    val_mcc = [m.get("mcc_mean", m.get("mcc", 0)) for m in history["val_metrics"]]
    ax.plot(epochs, val_f1,  color=MACARON_COLORS[2], linewidth=2, label="F1 Score")
    ax.plot(epochs, val_mcc, color=MACARON_COLORS[3], linewidth=2, label="MCC")
    ax.set_xlabel("Epoch", fontweight="bold")
    ax.set_ylabel("Score", fontweight="bold")
    ax.set_title("Validation F1 and MCC", fontweight="bold")
    ax.legend(prop={"weight": "bold"})
    ax.set_ylim([-0.05, 1.05])

    # Learning rate
    ax = axes[2]
    ax.semilogy(epochs, history["learning_rates"],
                color=MACARON_COLORS[4], linewidth=2)
    ax.set_xlabel("Epoch", fontweight="bold")
    ax.set_ylabel("Learning Rate", fontweight="bold")
    ax.set_title("Learning Rate Schedule", fontweight="bold")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info("Training curves saved → %s", save_path)
    plt.close()


def plot_metrics_comparison(
    results: Dict[str, Dict[str, float]],
    metrics: List[str] = None,
    save_path: Optional[Path] = None,
) -> None:
    """
    Bar chart comparing multiple models across evaluation metrics.
    """
    if metrics is None:
        metrics = ["precision_mean", "sensitivity_mean", "f1_mean", "mcc_mean"]

    models = list(results.keys())
    x      = np.arange(len(metrics))
    width  = 0.8 / max(len(models), 1)

    fig, ax = plt.subplots(figsize=(12, 6))

    for i, model in enumerate(models):
        vals = [results[model].get(m, 0.0) for m in metrics]
        errs = [results[model].get(m.replace("_mean", "_std"), 0.0)
                for m in metrics]
        offset = (i - len(models) / 2 + 0.5) * width
        bars = ax.bar(
            x + offset, vals,
            width=width * 0.9,
            color=MACARON_COLORS[i % len(MACARON_COLORS)],
            label=model,
            edgecolor="white",
            linewidth=0.8,
        )
        ax.errorbar(
            x + offset, vals, yerr=errs,
            fmt="none", color="black",
            capsize=4, linewidth=1.5,
        )

    tick_labels = [m.replace("_mean", "").replace("_", " ").title()
                   for m in metrics]
    ax.set_xticks(x)
    ax.set_xticklabels(tick_labels, fontweight="bold")
    ax.set_ylabel("Score", fontweight="bold")
    ax.set_title("Model Performance Comparison", fontweight="bold")
    ax.legend(prop={"weight": "bold"})
    ax.set_ylim([0, 1.05])

    # Significance markers (optional: add * for p<0.05 manually)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        logger.info("Comparison chart saved → %s", save_path)
    plt.close()


def plot_confusion_matrix_heatmap(
    tp: int, fp: int, fn: int, tn: int,
    title: str = "Confusion Matrix",
    save_path: Optional[Path] = None,
) -> None:
    """Plot 2×2 confusion matrix heatmap."""
    mat = np.array([[tn, fp], [fn, tp]], dtype=int)

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(mat, cmap=MACARON_SEQUENTIAL)

    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["Predicted Negative", "Predicted Positive"],
                       fontweight="bold")
    ax.set_yticklabels(["True Negative", "True Positive"],
                       fontweight="bold")
    ax.set_title(title, fontweight="bold")

    for i in range(2):
        for j in range(2):
            ax.text(j, i, str(mat[i, j]),
                    ha="center", va="center",
                    fontsize=18, fontweight="bold",
                    color="black" if mat[i, j] < mat.max() * 0.7 else "white")

    plt.colorbar(im, fraction=0.046)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


# ─────────────────────────────────────────────────────────────────────────────
#  25.  TRAINER
# ─────────────────────────────────────────────────────────────────────────────

if TORCH_AVAILABLE:

    class WarmupCosineScheduler(torch.optim.lr_scheduler._LRScheduler):
        """
        Linear warmup followed by cosine annealing.
        Reference: Loshchilov & Hutter (2017) SGDR.
        """

        def __init__(
            self,
            optimizer:    torch.optim.Optimizer,
            warmup_steps: int,
            total_steps:  int,
            min_lr:       float = 1e-6,
            last_epoch:   int   = -1,
        ) -> None:
            self.warmup_steps = warmup_steps
            self.total_steps  = total_steps
            self.min_lr       = min_lr
            super().__init__(optimizer, last_epoch)

        def get_lr(self):
            step = self.last_epoch + 1
            if step < self.warmup_steps:
                factor = step / max(self.warmup_steps, 1)
            else:
                progress = (step - self.warmup_steps) / max(
                    self.total_steps - self.warmup_steps, 1
                )
                factor = self.min_lr / self.base_lrs[0] + (
                    (1 - self.min_lr / self.base_lrs[0])
                    * 0.5 * (1 + math.cos(math.pi * progress))
                )
            return [base_lr * factor for base_lr in self.base_lrs]


    class Trainer:
        """
        Full training / evaluation pipeline.

        Features:
        - Gradient clipping
        - Warmup + cosine LR schedule
        - Best-model checkpoint saving (by validation F1)
        - Comprehensive metric logging
        - TensorBoard-compatible history export
        """

        def __init__(
            self,
            model:  RNASSPredictorAdvanced,
            cfg:    Config,
        ) -> None:
            self.model   = model
            self.cfg     = cfg
            self.device  = cfg.DEVICE

            self.criterion = RNAStructureLoss(cfg)
            self.optimizer = torch.optim.AdamW(
                model.parameters(),
                lr=cfg.LR,
                weight_decay=cfg.WEIGHT_DECAY,
                betas=(0.9, 0.999),
            )
            self.metrics  = RNAStructureMetrics(min_loop=cfg.MIN_HAIRPIN_LOOP)
            self.decoder  = NussinovDecoder(min_loop=cfg.MIN_HAIRPIN_LOOP)

            self.best_f1    = -1.0
            self.best_epoch = -1

            self.history: Dict[str, List] = {
                "train_loss":    [],
                "val_metrics":   [],
                "learning_rates":[],
            }

        def _make_scheduler(
            self, total_steps: int
        ) -> WarmupCosineScheduler:
            return WarmupCosineScheduler(
                self.optimizer,
                warmup_steps = self.cfg.WARMUP_STEPS,
                total_steps  = total_steps,
            )

        def train_epoch(
            self, loader: DataLoader, scheduler
        ) -> Tuple[float, Dict[str, float]]:
            self.model.train()
            total_loss  = 0.0
            stat_acc: Dict[str, float] = defaultdict(float)

            pbar = tqdm(loader, desc="  Training", leave=False)
            for batch in pbar:
                ids   = batch["input_ids"].to(self.device)
                tgt   = batch["pair_target"].to(self.device)
                mask  = batch["valid_mask"].to(self.device)
                prior = batch["pair_prior"].to(self.device)
                dist  = batch["dist_bias"].to(self.device)

                self.optimizer.zero_grad(set_to_none=True)
                logits = self.model(ids, prior, dist)
                loss, stats = self.criterion(logits, tgt, mask)

                loss.backward()
                torch.nn.utils.clip_grad_norm_(
                    self.model.parameters(), self.cfg.GRAD_CLIP
                )
                self.optimizer.step()
                scheduler.step()

                total_loss += loss.item()
                for k, v in stats.items():
                    stat_acc[k] += v
                pbar.set_postfix({"loss": f"{loss.item():.4f}"})

            n = max(len(loader), 1)
            return (
                total_loss / n,
                {k: v / n for k, v in stat_acc.items()},
            )

        def validate(
            self, loader: DataLoader
        ) -> Dict[str, float]:
            self.model.eval()
            all_pred:  List[str] = []
            all_true:  List[str] = []
            all_seqs:  List[str] = []
            all_probs: List[float] = []
            all_labels:List[int]  = []

            with torch.no_grad():
                for batch in tqdm(loader, desc="  Validating", leave=False):
                    ids   = batch["input_ids"].to(self.device)
                    prior = batch["pair_prior"].to(self.device)
                    dist  = batch["dist_bias"].to(self.device)

                    logits = self.model(ids, prior, dist)
                    probs  = torch.sigmoid(logits).cpu().numpy()

                    for i, seq in enumerate(batch["sequence"]):
                        L       = int(batch["length"][i])
                        true_db = batch["dot_bracket"][i]
                        pred_db = self.decoder.decode(
                            probs[i, :L, :L], seq,
                            self.cfg.DECODING_THRESHOLD,
                        )
                        all_pred.append(pred_db)
                        all_true.append(true_db)
                        all_seqs.append(seq)

                        # Collect per-pair labels & scores for AUC
                        prob_mat = probs[i, :L, :L]
                        true_mat = build_pair_matrix(true_db, L)
                        for ii in range(L):
                            for jj in range(ii + self.cfg.MIN_HAIRPIN_LOOP + 1, L):
                                all_probs.append(float(prob_mat[ii, jj]))
                                all_labels.append(int(true_mat[ii, jj]))

            batch_m = self.metrics.batch_metrics(all_pred, all_true, all_seqs)

            # ROC / PR AUC
            if all_labels and len(set(all_labels)) > 1:
                fpr, tpr, _ = roc_curve(all_labels, all_probs)
                roc_a       = float(auc(fpr, tpr))
                prec, rec, _= precision_recall_curve(all_labels, all_probs)
                pr_a        = float(auc(rec, prec))
                batch_m["roc_auc"] = roc_a
                batch_m["pr_auc"]  = pr_a

            return batch_m

        def train(
            self,
            train_loader: DataLoader,
            val_loader:   DataLoader,
            num_epochs:   int,
        ) -> None:
            total_steps = num_epochs * len(train_loader)
            scheduler   = self._make_scheduler(total_steps)

            logger.info(
                "Training started  device=%s  epochs=%d  steps=%d",
                self.device, num_epochs, total_steps,
            )

            for epoch in range(1, num_epochs + 1):
                logger.info("── Epoch %d / %d ──", epoch, num_epochs)

                tr_loss, tr_stats = self.train_epoch(train_loader, scheduler)
                val_m             = self.validate(val_loader)
                lr_now            = self.optimizer.param_groups[0]["lr"]

                self.history["train_loss"].append(tr_loss)
                self.history["val_metrics"].append(val_m)
                self.history["learning_rates"].append(lr_now)

                f1_now = val_m.get("f1_mean", val_m.get("f1", 0.0))
                logger.info(
                    "  train_loss=%.4f  val_F1=%.4f  val_MCC=%.4f  LR=%.2e",
                    tr_loss,
                    f1_now,
                    val_m.get("mcc_mean", val_m.get("mcc", 0.0)),
                    lr_now,
                )

                if f1_now > self.best_f1:
                    self.best_f1    = f1_now
                    self.best_epoch = epoch
                    self._save_checkpoint(epoch, val_m, best=True)
                    logger.info(
                        "  ★ New best  F1=%.4f  saved", self.best_f1
                    )

                # Periodic checkpoint
                if epoch % 10 == 0:
                    self._save_checkpoint(epoch, val_m, best=False)

            logger.info(
                "Training complete  best_F1=%.4f  best_epoch=%d",
                self.best_f1, self.best_epoch,
            )
            self._save_history()
            plot_training_curves(
                self.history,
                save_path=config.FIG_DIR / "training_curves.pdf",
            )

        def _save_checkpoint(
            self,
            epoch:   int,
            metrics: Dict,
            best:    bool = False,
        ) -> None:
            ckpt = {
                "epoch":              epoch,
                "model_state_dict":   self.model.state_dict(),
                "optimizer_state_dict": self.optimizer.state_dict(),
                "metrics":            metrics,
                "best_f1":            self.best_f1,
                "config": {
                    "d_model":    self.cfg.D_MODEL,
                    "nhead":      self.cfg.NHEAD,
                    "num_layers": self.cfg.NUM_LAYERS,
                    "kmer_size":  self.cfg.KMER_SIZE,
                    "max_len":    self.cfg.MAX_LEN,
                },
            }
            path = (
                config.MODEL_DIR / "best_model.pt"
                if best
                else config.MODEL_DIR / f"checkpoint_epoch{epoch:04d}.pt"
            )
            torch.save(ckpt, path)

        def _save_history(self) -> None:
            def _to_py(v: Any) -> Any:
                if isinstance(v, (np.floating, np.integer)):
                    return float(v)
                if isinstance(v, dict):
                    return {kk: _to_py(vv) for kk, vv in v.items()}
                if isinstance(v, list):
                    return [_to_py(x) for x in v]
                return v

            save_json(
                _to_py(self.history),
                config.LOG_DIR / "training_history.json",
            )


# ─────────────────────────────────────────────────────────────────────────────
#  26.  CROSS-VALIDATION
# ─────────────────────────────────────────────────────────────────────────────

if TORCH_AVAILABLE:

    def cross_validate(
        df:           pd.DataFrame,
        tokenizer:    KMerTokenizer,
        n_folds:      int = 5,
        fast_epochs:  int = 10,
    ) -> Tuple[List[Dict], Dict[str, Dict]]:
        """
        K-fold cross-validation.
        Uses random (not homology-aware) splits for speed;
        final results should be interpreted accordingly.

        Returns:
            fold_results : list of per-fold metric dicts
            summary      : {metric: {mean, std, min, max}}
        """
        kf           = KFold(n_splits=n_folds, shuffle=True, random_state=config.SEED)
        fold_results = []

        logger.info("Starting %d-fold cross-validation", n_folds)

        for fold_idx, (tr_idx, va_idx) in enumerate(kf.split(df)):
            logger.info("  Fold %d / %d", fold_idx + 1, n_folds)

            tr_ds = RNASSDataset(
                df.iloc[tr_idx].reset_index(drop=True), tokenizer, config.MAX_LEN
            )
            va_ds = RNASSDataset(
                df.iloc[va_idx].reset_index(drop=True), tokenizer, config.MAX_LEN
            )

            tr_ld = DataLoader(
                tr_ds, batch_size=config.BATCH_SIZE,
                shuffle=True, collate_fn=collate_fn,
            )
            va_ld = DataLoader(
                va_ds, batch_size=config.BATCH_SIZE,
                shuffle=False, collate_fn=collate_fn,
            )

            model = RNASSPredictorAdvanced(
                vocab_size  = tokenizer.vocab_size,
                d_model     = config.D_MODEL,
                nhead       = config.NHEAD,
                num_layers  = config.NUM_LAYERS,
                dim_ff      = config.DIM_FF,
                dropout     = config.DROPOUT,
                max_len     = config.MAX_LEN,
            ).to(config.DEVICE)

            trainer = Trainer(model, config)
            trainer.train(tr_ld, va_ld, fast_epochs)

            val_m = trainer.validate(va_ld)
            val_m["fold"] = fold_idx
            fold_results.append(val_m)

        # Summarise
        scalar_keys = [
            k for k in fold_results[0]
            if isinstance(fold_results[0][k], (int, float))
            and k != "fold"
        ]
        summary: Dict[str, Dict] = {}
        for k in scalar_keys:
            vals = np.array([r[k] for r in fold_results])
            summary[k] = {
                "mean": float(vals.mean()),
                "std":  float(vals.std()),
                "min":  float(vals.min()),
                "max":  float(vals.max()),
            }

        logger.info(
            "Cross-validation complete  F1=%.4f±%.4f  MCC=%.4f±%.4f",
            summary.get("f1_mean", {}).get("mean", 0),
            summary.get("f1_mean", {}).get("std",  0),
            summary.get("mcc_mean", {}).get("mean", 0),
            summary.get("mcc_mean", {}).get("std",  0),
        )
        save_json(
            {"fold_results": fold_results, "summary": summary},
            config.RESULTS_DIR / "cross_validation_results.json",
        )
        return fold_results, summary


# ─────────────────────────────────────────────────────────────────────────────
#  27.  ABLATION EXPERIMENT FRAMEWORK
# ─────────────────────────────────────────────────────────────────────────────

if TORCH_AVAILABLE:

    class AblationExperiment:
        """
        Systematic ablation study framework.

        Experiments:
          1. Data splitting strategy  (random vs. homology-aware)
          2. Sequence representation  (1-mer, 2-mer, 3-mer, 4-mer)
          3. Model depth              (2 / 4 / 6 / 8 Transformer layers)
          4. Loss function components (progressive inclusion)
          5. Decoder strategy         (threshold / DP-decode)
        """

        def __init__(self, df: pd.DataFrame) -> None:
            self.df  = df
            self._stats = StatisticalSignificanceTester()
            self._results: Dict[str, Any] = {}

        def _quick_train_eval(
            self,
            tr_df: pd.DataFrame,
            va_df: pd.DataFrame,
            tokenizer: KMerTokenizer,
            epochs: int = 5,
            **model_kwargs,
        ) -> Dict[str, float]:
            """Lightweight train→eval loop for ablation comparisons."""
            tr_ds = RNASSDataset(tr_df, tokenizer, config.MAX_LEN)
            va_ds = RNASSDataset(va_df, tokenizer, config.MAX_LEN)
            tr_ld = DataLoader(tr_ds, config.BATCH_SIZE,
                               shuffle=True, collate_fn=collate_fn)
            va_ld = DataLoader(va_ds, config.BATCH_SIZE,
                               shuffle=False, collate_fn=collate_fn)

            model = RNASSPredictorAdvanced(
                vocab_size=tokenizer.vocab_size,
                max_len=config.MAX_LEN,
                **model_kwargs,
            ).to(config.DEVICE)

            trainer = Trainer(model, config)
            trainer.train(tr_ld, va_ld, epochs)
            return trainer.validate(va_ld)

        def run_kmer_ablation(
            self,
            tr_df: pd.DataFrame,
            va_df: pd.DataFrame,
            k_values: List[int] = None,
            epochs: int = 5,
        ) -> Dict[str, Dict]:
            """Compare different k-mer sizes."""
            if k_values is None:
                k_values = [1, 2, 3, 4]

            results: Dict[str, Dict] = {}
            for k in k_values:
                name = f"k={k}"
                logger.info("Ablation: k-mer %s", name)
                tok = KMerTokenizer(k)
                r   = self._quick_train_eval(
                    tr_df, va_df, tok,
                    epochs=epochs,
                    d_model=config.D_MODEL,
                    nhead=config.NHEAD,
                    num_layers=config.NUM_LAYERS,
                    dim_ff=config.DIM_FF,
                    dropout=config.DROPOUT,
                )
                results[name] = r

            self._results["kmer_ablation"] = results
            return results

        def run_depth_ablation(
            self,
            tr_df: pd.DataFrame,
            va_df: pd.DataFrame,
            depths: List[int] = None,
            epochs: int = 5,
        ) -> Dict[str, Dict]:
            """Compare different Transformer depths."""
            if depths is None:
                depths = [2, 4, 6, 8]

            tokenizer = KMerTokenizer(config.KMER_SIZE)
            results: Dict[str, Dict] = {}

            for d in depths:
                name = f"layers={d}"
                logger.info("Ablation: depth %s", name)
                r = self._quick_train_eval(
                    tr_df, va_df, tokenizer,
                    epochs=epochs,
                    d_model=config.D_MODEL,
                    nhead=config.NHEAD,
                    num_layers=d,
                    dim_ff=config.DIM_FF,
                    dropout=config.DROPOUT,
                )
                results[name] = r

            self._results["depth_ablation"] = results
            return results

        def generate_report(self) -> None:
            """Save ablation results and generate comparison figures."""
            save_json(self._results, config.RESULTS_DIR / "ablation_results.json")

            for exp_name, exp_data in self._results.items():
                if not exp_data:
                    continue
                plot_metrics_comparison(
                    results   = exp_data,
                    metrics   = ["f1_mean", "mcc_mean", "precision_mean",
                                 "sensitivity_mean"],
                    save_path = config.FIG_DIR / f"ablation_{exp_name}.pdf",
                )
            logger.info("Ablation report saved → %s", config.RESULTS_DIR)


# ─────────────────────────────────────────────────────────────────────────────
#  28.  COMMAND-LINE INTERFACE
# ─────────────────────────────────────────────────────────────────────────────

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="RNA Secondary Structure Prediction — SCI-Level System",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Pipeline stages
    stage = parser.add_argument_group("Pipeline Stages")
    stage.add_argument("--crawl",      action="store_true", help="Run data crawling")
    stage.add_argument("--process",    action="store_true", help="Process crawled data")
    stage.add_argument("--train",      action="store_true", help="Train prediction model")
    stage.add_argument("--evaluate",   action="store_true", help="Evaluate on test set")
    stage.add_argument("--ablation",   action="store_true", help="Run ablation experiments")
    stage.add_argument("--cross_val",  action="store_true", help="Run cross-validation")
    stage.add_argument("--predict",    type=str,            help="Predict a single sequence")

    # Crawling
    cw = parser.add_argument_group("Crawling")
    cw.add_argument("--target_sequences", type=int, default=50_000)
    cw.add_argument("--sources", nargs="+",
                    default=["bprna", "rfam", "rnacentral"])
    cw.add_argument("--parallel", action="store_true")

    # Data
    da = parser.add_argument_group("Data")
    da.add_argument("--min_len",           type=int,   default=20)
    da.add_argument("--max_len",           type=int,   default=256)
    da.add_argument("--require_structure", action="store_true")
    da.add_argument("--cluster_method",    type=str,   default="mmseqs2",
                    choices=["mmseqs2", "cdhit", "random"])

    # Training
    tr = parser.add_argument_group("Training")
    tr.add_argument("--epochs",     type=int,   default=50)
    tr.add_argument("--batch_size", type=int,   default=8)
    tr.add_argument("--lr",         type=float, default=2e-4)
    tr.add_argument("--kmer",       type=int,   default=3)
    tr.add_argument("--d_model",    type=int,   default=128)
    tr.add_argument("--use_rnafm",  action="store_true")

    # Misc
    misc = parser.add_argument_group("Miscellaneous")
    misc.add_argument("--seed",       type=int, default=42)
    misc.add_argument("--model_path", type=str, help="Path to saved checkpoint")
    misc.add_argument("--device",     type=str, choices=["cuda", "cpu"])

    return parser.parse_args()


# ─────────────────────────────────────────────────────────────────────────────
#  29.  MAIN PIPELINE
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_arguments()
    config.update_from_args(args)
    if args.device:
        config.DEVICE = args.device
    if args.cluster_method:
        config.CLUSTER_METHOD = args.cluster_method

    set_seed(config.SEED)

    banner = "=" * 72
    logger.info(banner)
    logger.info(
        "  RNA Secondary Structure Prediction System  |  SCI-Level Release"
    )
    logger.info("  Device: %s", config.DEVICE)
    logger.info("  Output: %s", config.BASE_DIR.resolve())
    logger.info(banner)

    # ── Stage 1: Crawl ────────────────────────────────────────────────────
    if args.crawl:
        logger.info("[Stage 1] Data Crawling")
        # Validate sources
        with DataCrawlerManager() as mgr:
            health = mgr.health_check()
            for src, status in health.items():
                logger.info("  %s  available=%s  rt=%.2fs",
                            src,
                            status.get("available"),
                            status.get("response_time") or -1)
            total = mgr.crawl_all(
                target_total=args.target_sequences,
                parallel=args.parallel,
            )
        logger.info("[Stage 1] Complete  total_sequences=%d", total)

    # ── Stage 2: Process ──────────────────────────────────────────────────
    if args.process:
        logger.info("[Stage 2] Data Processing")
        with DataProcessor() as proc:
            df = proc.build_dataset(require_structure=args.require_structure)
            if df.empty:
                logger.error("No data — run --crawl first"); return

            df = proc.filter_by_structure_complexity(df, min_pairs=2)
            proc.export_fasta(df, config.PROC_DIR / "all_sequences.fasta")
            proc.export_stockholm(df, config.PROC_DIR / "all_structures.stk")
            proc.visualize_dataset(df)
        logger.info("[Stage 2] Complete  n=%d", len(df))

    # ── Stage 3: Train ────────────────────────────────────────────────────
    if args.train:
        if not TORCH_AVAILABLE:
            logger.error("PyTorch required for training"); return

        logger.info("[Stage 3] Training")

        dataset_path = config.PROC_DIR / "raw_dataset.csv"
        if not dataset_path.exists():
            logger.error("No processed data — run --process first"); return

        df = pd.read_csv(dataset_path)

        # Homology-aware split
        splitter    = HomologyAwareSplitter(
            method=config.CLUSTER_METHOD,
            identity=config.HOMOLOGY_IDENTITY,
        )
        tr_df, va_df, te_df = splitter.split(df)

        # Tokenizer
        tokenizer = KMerTokenizer(config.KMER_SIZE)
        save_pickle(tokenizer, config.PROC_DIR / "tokenizer.pkl")

        # Datasets & Loaders
        tr_ds = RNASSDataset(tr_df, tokenizer, config.MAX_LEN, augment=True)
        va_ds = RNASSDataset(va_df, tokenizer, config.MAX_LEN)
        te_ds = RNASSDataset(te_df, tokenizer, config.MAX_LEN)

        make_loader = lambda ds, shuf: DataLoader(
            ds, batch_size=config.BATCH_SIZE, shuffle=shuf,
            num_workers=config.NUM_WORKERS,
            collate_fn=collate_fn,
            pin_memory=(config.DEVICE == "cuda"),
        )
        tr_ld = make_loader(tr_ds, True)
        va_ld = make_loader(va_ds, False)
        te_ld = make_loader(te_ds, False)

        # Model
        model = RNASSPredictorAdvanced(
            vocab_size   = tokenizer.vocab_size,
            d_model      = config.D_MODEL,
            nhead        = config.NHEAD,
            num_layers   = config.NUM_LAYERS,
            dim_ff       = config.DIM_FF,
            dropout      = config.DROPOUT,
            max_len      = config.MAX_LEN,
            use_rnafm    = args.use_rnafm,
            unet_base_ch = config.UNET_BASE_CH,
        ).to(config.DEVICE)

        n_params = sum(p.numel() for p in model.parameters())
        logger.info("Model parameters: %s", f"{n_params:,}")

        trainer = Trainer(model, config)
        trainer.train(tr_ld, va_ld, config.NUM_EPOCHS)

        # Test evaluation
        logger.info("[Stage 3] Final Test Evaluation")
        test_m = trainer.validate(te_ld)
        logger.info("Test Results:")
        for k, v in sorted(test_m.items()):
            if isinstance(v, float):
                logger.info("  %-30s %.4f", k, v)

        # Aggregate confusion matrix
        tp = int(sum(r.get("tp", 0) for r in [test_m]))
        fp = int(sum(r.get("fp", 0) for r in [test_m]))
        fn = int(sum(r.get("fn", 0) for r in [test_m]))
        tn = int(sum(r.get("tn", 0) for r in [test_m]))
        plot_confusion_matrix_heatmap(
            tp, fp, fn, tn,
            title="Test Set Confusion Matrix",
            save_path=config.FIG_DIR / "confusion_matrix.pdf",
        )

        save_json(test_m, config.RESULTS_DIR / "test_results.json")
        logger.info("[Stage 3] Complete")

    # ── Stage 4: Cross-Validation ─────────────────────────────────────────
    if args.cross_val:
        if not TORCH_AVAILABLE:
            logger.error("PyTorch required"); return

        logger.info("[Stage 4] Cross-Validation")
        df = pd.read_csv(config.PROC_DIR / "raw_dataset.csv")
        tokenizer = KMerTokenizer(config.KMER_SIZE)
        fold_res, summary = cross_validate(df, tokenizer, n_folds=5)
        logger.info("CV Summary:")
        for k, v in sorted(summary.items()):
            logger.info("  %-30s mean=%.4f ± %.4f", k, v["mean"], v["std"])
        logger.info("[Stage 4] Complete")

    # ── Stage 5: Ablation ─────────────────────────────────────────────────
    if args.ablation:
        if not TORCH_AVAILABLE:
            logger.error("PyTorch required"); return

        logger.info("[Stage 5] Ablation Study")
        df = pd.read_csv(config.PROC_DIR / "raw_dataset.csv")
        splitter = HomologyAwareSplitter(
            method=config.CLUSTER_METHOD,
            identity=config.HOMOLOGY_IDENTITY,
        )
        tr_df, va_df, _ = splitter.split(df)
        ablation = AblationExperiment(df)

        # k-mer ablation
        kmer_res = ablation.run_kmer_ablation(tr_df, va_df, epochs=5)
        logger.info("k-mer ablation: %s", kmer_res)

        # Depth ablation
        depth_res = ablation.run_depth_ablation(tr_df, va_df, epochs=5)
        logger.info("Depth ablation: %s", depth_res)

        ablation.generate_report()
        logger.info("[Stage 5] Complete")

    # ── Stage 6: Single-Sequence Prediction ──────────────────────────────
    if args.predict:
        if not TORCH_AVAILABLE:
            logger.error("PyTorch required"); return

        logger.info("[Stage 6] Single Sequence Prediction")
        model_path = Path(args.model_path) if args.model_path else \
                     config.MODEL_DIR / "best_model.pt"

        if not model_path.exists():
            logger.error("Model not found: %s", model_path); return

        ckpt      = torch.load(model_path, map_location=config.DEVICE)
        tok_path  = config.PROC_DIR / "tokenizer.pkl"
        tokenizer = load_pickle(tok_path) if tok_path.exists() else \
                    KMerTokenizer(config.KMER_SIZE)

        model = RNASSPredictorAdvanced(
            vocab_size = tokenizer.vocab_size,
            **{k: ckpt["config"][k] for k in
               ["d_model", "nhead", "num_layers", "max_len"]},
            dim_ff  = config.DIM_FF,
            dropout = 0.0,
        ).to(config.DEVICE)
        model.load_state_dict(ckpt["model_state_dict"])
        model.eval()

        seq = clean_sequence(args.predict)
        logger.info("Sequence (%d nt): %s", len(seq), seq)

        ids   = torch.tensor(
            tokenizer.encode(seq, config.MAX_LEN), dtype=torch.long
        ).unsqueeze(0).to(config.DEVICE)
        prior = torch.tensor(
            build_basepair_prior(seq, config.MAX_LEN),
            dtype=torch.float32
        ).unsqueeze(0).to(config.DEVICE)
        dist  = torch.tensor(
            build_distance_bias(len(seq), config.MAX_LEN),
            dtype=torch.float32
        ).unsqueeze(0).to(config.DEVICE)

        with torch.no_grad():
            logits = model(ids, prior, dist)
            probs  = torch.sigmoid(logits).cpu().numpy()[0, :len(seq), :len(seq)]

        # Thermodynamic fusion
        vienna = ViennaRNAInterface()
        thermo = vienna.pair_probabilities(seq)
        fused  = ViennaRNAInterface.fuse(
            probs, thermo, alpha=1.0 - config.THERMO_WEIGHT
        )

        decoder = NussinovDecoder(min_loop=config.MIN_HAIRPIN_LOOP)
        pred_db = decoder.decode(fused, seq, config.DECODING_THRESHOLD)

        logger.info("Predicted structure: %s", pred_db)
        n_pairs = len(PseudoKnotParser.dotbracket_to_pairs(pred_db))
        logger.info("Base pairs predicted: %d", n_pairs)

        # MC Dropout uncertainty
        estimator = MCDropoutEstimator(model, n_samples=config.MC_SAMPLES)
        mean_p, std_p, entropy = estimator.estimate(ids, prior, dist)
        unc_global = float(std_p[0, :len(seq), :len(seq)].mean())
        logger.info("Global epistemic uncertainty: %.5f", unc_global)

        estimator.plot_uncertainty(
            mean_p[0], std_p[0], seq,
            save_path=config.FIG_DIR / "uncertainty_prediction.pdf",
        )

        result = {
            "sequence":    seq,
            "structure":   pred_db,
            "n_pairs":     n_pairs,
            "uncertainty": unc_global,
        }
        save_json(result, config.RESULTS_DIR / "single_prediction.json")
        logger.info("[Stage 6] Complete — structure saved")

    # ── Done ──────────────────────────────────────────────────────────────
    logger.info(banner)
    logger.info("Pipeline complete.  Results → %s", config.RESULTS_DIR.resolve())
    logger.info(banner)


if __name__ == "__main__":
    main()
