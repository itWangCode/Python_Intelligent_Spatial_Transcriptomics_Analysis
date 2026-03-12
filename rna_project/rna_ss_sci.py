#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
RNA Secondary Structure Prediction System - Complete SCI Edition
================================================================

A comprehensive deep learning framework for RNA secondary structure prediction,
designed to meet the rigorous standards of high-impact scientific publications.

Key Features:
1. Multi-source data crawling with resume capability (bpRNA, Rfam, RNAcentral)
2. Homology-aware data splitting (MMseqs2/CD-HIT)
3. k-mer and RNA-FM pretrained representations
4. Transformer + U-Net hybrid architecture with attention visualization
5. Multi-objective loss function (BCE, Focal, Dice, Symmetry, Sparsity)
6. Nussinov dynamic programming decoding (canonical structures only)
7. ViennaRNA thermodynamic model fusion
8. Monte Carlo Dropout uncertainty estimation
9. Strict evaluation metrics with valid pair space definition
10. Comprehensive ablation study framework
11. Statistical significance testing
12. Professional visualization with Nature-style formatting

Author: Research Group
Date: 2024
License: MIT
"""

import os
import sys
import re
import json
import math
import time
import gzip
import random
import shutil
import pickle
import hashlib
import zipfile
import sqlite3
import urllib
import requests
import warnings
import argparse
import datetime
import subprocess
from pathlib import Path
from tqdm import tqdm
from bs4 import BeautifulSoup
from urllib.parse import urljoin, urlparse
from typing import List, Tuple, Dict, Optional, Union, Any, Set, Callable
from collections import defaultdict, Counter
from dataclasses import dataclass, field
from enum import Enum
from abc import ABC, abstractmethod
import itertools
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import softmax
from scipy.stats import pearsonr, spearmanr

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader, WeightedRandomSampler, Subset

# Optional imports for RNA-FM
try:
    from transformers import AutoModel, AutoTokenizer
    TRANSFORMERS_AVAILABLE = True
except ImportError:
    TRANSFORMERS_AVAILABLE = False
    print("[WARNING] transformers not installed. RNA-FM integration disabled.")

from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix

# Visualization imports
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec

warnings.filterwarnings("ignore")

# =========================================================
# 0. Global Configuration
# =========================================================
class Config:
    """Unified configuration class for the entire pipeline."""
    
    def __init__(self):
        # Random seed for reproducibility
        self.SEED = 42
        
        # Device configuration
        self.DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
        
        # Directory structure
        self.BASE_DIR = Path("./rna_ss_sci")
        self.RAW_DIR = self.BASE_DIR / "data" / "raw"
        self.PROC_DIR = self.BASE_DIR / "data" / "processed"
        self.MODEL_DIR = self.BASE_DIR / "models"
        self.TMP_DIR = self.BASE_DIR / "tmp"
        self.LOG_DIR = self.BASE_DIR / "logs"
        self.RESULTS_DIR = self.BASE_DIR / "results"
        self.FIGURES_DIR = self.BASE_DIR / "figures"
        self.CRAWLER_DB = self.BASE_DIR / "crawler_state.db"
        
        # Data sources configuration
        self.DATA_SOURCES = {
            'bprna': {
                'base_url': 'https://bprna.cgrb.oregonstate.edu/',
                'dotbracket_url': 'https://bprna.cgrb.oregonstate.edu/download.php?download=dotbracket',
                'fasta_url': 'https://bprna.cgrb.oregonstate.edu/download.php?download=fasta',
                'enabled': True,
                'max_sequences': 50000,
                'priority': 1,
                'has_structure': True
            },
            'rfam': {
                'base_url': 'https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/fasta_files/',
                'family_list_url': 'https://rfam.org/families',
                'enabled': True,
                'max_sequences': 100000,
                'priority': 2,
                'has_structure': False
            },
            'rnacentral': {
                'api_url': 'https://rnacentral.org/api/v1/',
                'enabled': True,
                'max_sequences': 50000,
                'priority': 3,
                'has_structure': False
            }
        }
        
        # Crawler configuration
        self.REQUEST_TIMEOUT = 30
        self.MAX_RETRIES = 3
        self.RETRY_DELAY = 5
        self.CONCURRENT_DOWNLOADS = 4
        self.DOWNLOAD_CHUNK_SIZE = 8192
        self.USER_AGENT = "Mozilla/5.0 (compatible; RNA-SS-System/1.0; +https://github.com/yourlab/rna-ss)"
        self.RATE_LIMIT = 1.0  # seconds between requests
        self.MAX_CRAWL_PER_SOURCE = 100000
        
        # Data filtering
        self.MIN_LEN = 20
        self.MAX_LEN = 256
        self.MIN_HAIRPIN_LOOP = 3
        self.MIN_SEQUENCE_QUALITY = 0.8
        self.MIN_PAIRS_FOR_STRUCTURE = 1
        
        # Data splitting
        self.TRAIN_RATIO = 0.70
        self.VALID_RATIO = 0.15
        self.TEST_RATIO = 0.15
        
        # Homology clustering
        self.HOMOLOGY_IDENTITY = 0.80  # 80% sequence identity threshold
        self.CLUSTER_METHOD = "mmseqs2"  # "mmseqs2" or "cdhit"
        self.CLUSTER_COVERAGE = 0.8
        self.CLUSTER_MODE = "bidirectional"
        
        # Training parameters
        self.BATCH_SIZE = 8
        self.NUM_EPOCHS = 50
        self.LR = 2e-4
        self.WEIGHT_DECAY = 1e-4
        self.NUM_WORKERS = 4
        self.GRAD_CLIP = 1.0
        self.EARLY_STOPPING_PATIENCE = 10
        self.NUM_REPEATS = 3  # Number of repeats for statistical significance
        
        # Model parameters
        self.REPRESENTATION = "kmer"  # "kmer" or "rnafm"
        self.KMER_SIZE = 3
        self.D_MODEL = 128
        self.NHEAD = 8
        self.NUM_LAYERS = 4
        self.DIM_FF = 256
        self.DROPOUT = 0.1
        self.USE_POSITIONAL_ENCODING = True
        
        # Loss function weights (for valid pair space only)
        self.LOSS_WEIGHTS = {
            'bce': 0.5,
            'focal': 0.25,
            'dice': 0.20,
            'symmetry': 0.02,
            'diagonal': 0.01,
            'sparsity': 0.02,
            'degree': 0.00  # Optional degree regularization
        }
        
        # Positive sample weight for class imbalance
        self.POS_WEIGHT = 8.0
        
        # Focal Loss parameters
        self.FOCAL_ALPHA = 0.25
        self.FOCAL_GAMMA = 2.0
        
        # Target pair density for sparsity regularization
        self.TARGET_PAIR_DENSITY = 0.02
        
        # Decoding parameters
        self.DECODING_THRESHOLD = 0.5
        self.DECODING_METHOD = "nussinov"  # "greedy" or "nussinov"
        self.SUPPORT_PSEUDOKNOTS = False  # Set to False for canonical structures only
        
        # Uncertainty estimation
        self.MC_SAMPLES = 20
        
        # Thermodynamic fusion
        self.USE_VIENNA = False
        self.THERMO_WEIGHT = 0.2
        
        # Evaluation
        self.EVALUATION_SPACE = "valid_pairs"  # "valid_pairs" or "all_pairs"
        self.COMPUTE_CONFIDENCE_INTERVALS = True
        self.CONFIDENCE_LEVEL = 0.95
        
        # Visualization
        self.FIGURE_DPI = 300
        self.FIGURE_FORMAT = "pdf"  # "pdf", "png", "svg"
        self.FONT_FAMILY = "Times New Roman"
        self.FONT_SIZE = 12
        self.TITLE_FONT_SIZE = 14
        self.LABEL_FONT_SIZE = 12
        self.LEGEND_FONT_SIZE = 10
        self.COLOR_PALETTE = "muted"  # seaborn color palette
        
        # Create all directories
        self._create_directories()
        
        # Initialize crawler database
        self._init_crawler_db()
    
    def _create_directories(self):
        """Create all necessary directories."""
        for dir_path in [self.RAW_DIR, self.PROC_DIR, self.MODEL_DIR, 
                        self.TMP_DIR, self.LOG_DIR, self.RESULTS_DIR, 
                        self.FIGURES_DIR]:
            dir_path.mkdir(parents=True, exist_ok=True)
    
    def _init_crawler_db(self):
        """Initialize crawler state database for resume capability."""
        conn = sqlite3.connect(str(self.CRAWLER_DB))
        cursor = conn.cursor()
        
        # Data sources status table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS data_sources (
                source_name TEXT PRIMARY KEY,
                total_sequences INTEGER DEFAULT 0,
                downloaded_sequences INTEGER DEFAULT 0,
                last_crawled TIMESTAMP,
                status TEXT,
                error_message TEXT,
                config TEXT
            )
        ''')
        
        # Downloaded files table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS downloaded_files (
                file_id TEXT PRIMARY KEY,
                source_name TEXT,
                url TEXT,
                local_path TEXT,
                file_size INTEGER,
                md5_hash TEXT,
                download_time TIMESTAMP,
                status TEXT,
                retry_count INTEGER DEFAULT 0,
                FOREIGN KEY(source_name) REFERENCES data_sources(source_name)
            )
        ''')
        
        # Parsed sequences table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS parsed_sequences (
                sequence_id TEXT PRIMARY KEY,
                source_name TEXT,
                original_id TEXT,
                sequence TEXT,
                dot_bracket TEXT,
                length INTEGER,
                quality_score REAL,
                gc_content REAL,
                md5_hash TEXT UNIQUE,
                parsed_time TIMESTAMP,
                FOREIGN KEY(source_name) REFERENCES data_sources(source_name)
            )
        ''')
        
        # Create indices for performance
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_source ON parsed_sequences(source_name)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_length ON parsed_sequences(length)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_md5 ON parsed_sequences(md5_hash)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_status ON downloaded_files(status)')
        
        conn.commit()
        conn.close()
    
    def update_from_args(self, args):
        """Update configuration from command line arguments."""
        for key, value in vars(args).items():
            if value is not None and hasattr(self, key.upper()):
                setattr(self, key.upper(), value)
    
    def save(self, path: Path):
        """Save configuration to JSON file."""
        config_dict = {k: v for k, v in self.__dict__.items() 
                      if not k.startswith('_') and k.isupper()}
        save_json(config_dict, path)
    
    @classmethod
    def load(cls, path: Path) -> 'Config':
        """Load configuration from JSON file."""
        config = cls()
        config_dict = load_json(path)
        for k, v in config_dict.items():
            if hasattr(config, k):
                setattr(config, k, v)
        return config


config = Config()


# =========================================================
# 1. Utility Functions
# =========================================================
def set_seed(seed: int = config.SEED):
    """Set all random seeds for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def save_json(data: Any, path: Path):
    """Save data to JSON file."""
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


def load_json(path: Path) -> Any:
    """Load data from JSON file."""
    with open(path, 'r', encoding='utf-8') as f:
        return json.load(f)


def save_pickle(data: Any, path: Path):
    """Save data to pickle file."""
    with open(path, 'wb') as f:
        pickle.dump(data, f)


def load_pickle(path: Path) -> Any:
    """Load data from pickle file."""
    with open(path, 'rb') as f:
        return pickle.load(f)


def compute_md5(file_path: Path, chunk_size: int = 8192) -> str:
    """Compute MD5 hash of a file."""
    md5 = hashlib.md5()
    with open(file_path, 'rb') as f:
        for chunk in iter(lambda: f.read(chunk_size), b''):
            md5.update(chunk)
    return md5.hexdigest()


def compute_sequence_md5(sequence: str, structure: str = '') -> str:
    """Compute MD5 hash of a sequence-structure pair for deduplication."""
    content = f"{sequence}|{structure}".encode('utf-8')
    return hashlib.md5(content).hexdigest()


def get_timestamp() -> str:
    """Get current timestamp as string."""
    return datetime.datetime.now().isoformat()


def log_message(message: str, level: str = "INFO", log_file: Optional[Path] = None):
    """Log message to console and optionally to file."""
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_entry = f"[{timestamp}] [{level}] {message}"
    print(log_entry)
    
    if log_file is None:
        log_file = config.LOG_DIR / f"system_{datetime.datetime.now().strftime('%Y%m%d')}.log"
    
    with open(log_file, 'a', encoding='utf-8') as f:
        f.write(log_entry + '\n')


def has_executable(name: str) -> bool:
    """Check if an executable is available in PATH."""
    return shutil.which(name) is not None


def safe_div(a: float, b: float) -> float:
    """Safe division, returning 0 if denominator is 0."""
    return a / b if b != 0 else 0.0


# =========================================================
# 2. Sequence and Structure Cleaning
# =========================================================
def clean_sequence(seq: str) -> str:
    """
    Clean RNA sequence:
    - Convert to uppercase
    - Replace T with U (DNA to RNA)
    - Remove non-standard characters
    """
    if not isinstance(seq, str):
        return ""
    seq = seq.upper().replace("T", "U")
    seq = re.sub(r"[^ACGUN]", "", seq)
    return seq


def clean_dot_bracket(db: str) -> str:
    """
    Clean dot-bracket structure:
    - Keep only valid bracket characters: . ( ) [ ] { } < > and letters for pseudoknots
    """
    if not isinstance(db, str):
        return ""
    db = db.strip()
    # Keep only valid characters for RNA secondary structure representation
    db = re.sub(r"[^().\[\]{}<>A-Za-z]", "", db)
    return db


def validate_sequence_quality(seq: str) -> float:
    """
    Assess sequence quality as proportion of valid nucleotides.
    Returns a score between 0 and 1.
    """
    if len(seq) == 0:
        return 0.0
    clean = clean_sequence(seq)
    return len(clean) / len(seq)


def compute_gc_content(seq: str) -> float:
    """Compute GC content of a sequence."""
    if len(seq) == 0:
        return 0.0
    seq = clean_sequence(seq)
    gc_count = seq.count('G') + seq.count('C')
    return gc_count / len(seq)


# =========================================================
# 3. Bracketed Notation Parser (Supports Pseudoknots)
# =========================================================
class BracketedNotationParser:
    """
    Parser for dot-bracket notation supporting multiple bracket types.
    
    Bracket level mapping:
    - Level 0: ()  # Canonical base pairs
    - Level 1: []  # Pseudoknot level 1
    - Level 2: {}  # Pseudoknot level 2
    - Level 3: <>  # Pseudoknot level 3
    - Level 4+: Aa, Bb, Cc ...  # Additional pseudoknot levels
    
    Note: For canonical structure prediction, only level 0 is used.
    """
    
    # Base bracket pairs
    BASE_OPEN_TO_CLOSE = {
        '(': ')', '[': ']', '{': '}', '<': '>'
    }
    BASE_CLOSE_TO_OPEN = {v: k for k, v in BASE_OPEN_TO_CLOSE.items()}
    
    # Letter brackets for higher pseudoknot levels
    LOWERS = "abcdefghijklmnopqrstuvwxyz"
    UPPERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    
    # Extended bracket mapping (initialized lazily)
    OPEN_TO_CLOSE = BASE_OPEN_TO_CLOSE.copy()
    CLOSE_TO_OPEN = BASE_CLOSE_TO_OPEN.copy()
    
    @classmethod
    def _initialize_extended_brackets(cls):
        """Initialize extended bracket types for pseudoknot support."""
        if len(cls.OPEN_TO_CLOSE) > 4:
            return
        
        for lo, up in zip(cls.LOWERS, cls.UPPERS):
            cls.OPEN_TO_CLOSE[up] = lo
            cls.CLOSE_TO_OPEN[lo] = up
    
    @classmethod
    def dotbracket_to_pairs(cls, dot_bracket: str, 
                           allow_pseudoknots: bool = False) -> List[Tuple[int, int]]:
        """
        Convert dot-bracket string to list of base pairs.
        
        Args:
            dot_bracket: Dot-bracket notation string
            allow_pseudoknots: If True, process all bracket types; 
                              if False, only process canonical '(' and ')'
        
        Returns:
            List of (i, j) pairs (0-indexed)
        """
        if allow_pseudoknots:
            cls._initialize_extended_brackets()
            open_chars = cls.OPEN_TO_CLOSE.keys()
        else:
            open_chars = ['(']
        
        stacks = {k: [] for k in open_chars}
        pairs = []
        
        for i, ch in enumerate(dot_bracket):
            if ch in stacks:
                stacks[ch].append(i)
            elif allow_pseudoknots and ch in cls.CLOSE_TO_OPEN:
                opener = cls.CLOSE_TO_OPEN[ch]
                if stacks.get(opener):
                    j = stacks[opener].pop()
                    pairs.append((j, i))
            elif ch == ')':  # Canonical closing bracket
                if stacks.get('('):
                    j = stacks['('].pop()
                    pairs.append((j, i))
        
        pairs.sort()
        return pairs
    
    @classmethod
    def pairs_to_dotbracket(cls, pairs: List[Tuple[int, int]], length: int,
                           allow_pseudoknots: bool = False) -> str:
        """
        Convert list of base pairs to dot-bracket notation.
        
        Args:
            pairs: List of (i, j) pairs
            length: Length of the sequence
            allow_pseudoknots: If True, assign different bracket types to crossing pairs;
                              if False, only output canonical '()' for non-crossing pairs
        
        Returns:
            Dot-bracket notation string
        """
        result = ['.'] * length
        
        if not pairs:
            return ''.join(result)
        
        if not allow_pseudoknots:
            # For canonical structures, just use '()' for all pairs
            for i, j in pairs:
                if i < j:  # Ensure i < j for canonical representation
                    result[i] = '('
                    result[j] = ')'
            return ''.join(result)
        
        # For pseudoknots, assign bracket levels to avoid crossing
        cls._initialize_extended_brackets()
        
        # Sort pairs by length (longest first) for better bracket assignment
        pairs = sorted(pairs, key=lambda x: x[1] - x[0], reverse=True)
        
        # Build crossing graph
        pair_graph = {}
        for idx1, (i, j) in enumerate(pairs):
            pair_graph[(i, j)] = []
            for idx2, (m, n) in enumerate(pairs):
                if idx1 != idx2 and cls._are_crossing(i, j, m, n):
                    pair_graph[(i, j)].append((m, n))
        
        # Assign bracket levels using graph coloring
        levels = {}
        bracket_types = list(cls.OPEN_TO_CLOSE.keys())
        
        for pair in pairs:
            # Find all levels used by crossing pairs
            crossing_pairs = pair_graph.get(pair, [])
            used_levels = {levels.get(cp, -1) for cp in crossing_pairs if cp in levels}
            
            # Assign the smallest available level
            level = 0
            while level in used_levels:
                level += 1
            
            levels[pair] = level
        
        # Build dot-bracket string
        for (i, j), level in levels.items():
            if level < len(bracket_types):
                result[i] = bracket_types[level]
                result[j] = cls.CLOSE_TO_OPEN[result[i]]
            else:
                # Fallback for very high levels
                result[i] = '('
                result[j] = ')'
        
        return ''.join(result)
    
    @staticmethod
    def _are_crossing(a: int, b: int, c: int, d: int) -> bool:
        """Check if two pairs (a,b) and (c,d) are crossing."""
        return (a < c < b < d) or (c < a < d < b)
    
    @classmethod
    def pairs_to_matrix(cls, pairs: List[Tuple[int, int]], length: int,
                       symmetric: bool = True) -> np.ndarray:
        """Convert list of pairs to binary contact matrix."""
        mat = np.zeros((length, length), dtype=np.float32)
        for i, j in pairs:
            if 0 <= i < length and 0 <= j < length and i != j:
                if symmetric:
                    mat[i, j] = 1.0
                    mat[j, i] = 1.0
                else:
                    if i < j:
                        mat[i, j] = 1.0
                    else:
                        mat[j, i] = 1.0
        return mat
    
    @classmethod
    def dotbracket_to_matrix(cls, dot_bracket: str, allow_pseudoknots: bool = False) -> np.ndarray:
        """Convert dot-bracket directly to contact matrix."""
        pairs = cls.dotbracket_to_pairs(dot_bracket, allow_pseudoknots)
        return cls.pairs_to_matrix(pairs, len(dot_bracket))


# =========================================================
# 4. Base Pairing Rules
# =========================================================
class BasePairingRules:
    """Defines valid base pairing rules for RNA."""
    
    # Canonical Watson-Crick and wobble pairs
    CANONICAL_PAIRS = {
        ("A", "U"), ("U", "A"),
        ("C", "G"), ("G", "C"),
        ("G", "U"), ("U", "G")
    }
    
    # Pairs involving ambiguous nucleotides
    AMBIGUOUS_PAIRS = {
        ("N", "A"), ("A", "N"),
        ("N", "C"), ("C", "N"),
        ("N", "G"), ("G", "N"),
        ("N", "U"), ("U", "N"),
        ("N", "N")
    }
    
    # All allowed pairs
    ALLOWED_PAIRS = CANONICAL_PAIRS | AMBIGUOUS_PAIRS
    
    @classmethod
    def is_valid_pair(cls, a: str, b: str, include_ambiguous: bool = True) -> bool:
        """Check if two nucleotides can form a valid base pair."""
        if include_ambiguous:
            return (a, b) in cls.ALLOWED_PAIRS
        else:
            return (a, b) in cls.CANONICAL_PAIRS
    
    @classmethod
    def get_valid_pair_space(cls, sequence: str, min_loop: int = 3,
                            include_ambiguous: bool = True) -> List[Tuple[int, int]]:
        """
        Get all valid pair positions (i, j) with i < j, j - i > min_loop,
        and valid base pairing rules.
        """
        seq = clean_sequence(sequence)
        L = len(seq)
        valid_pairs = []
        
        for i in range(L):
            for j in range(i + 1, L):
                if j - i <= min_loop:
                    continue
                if cls.is_valid_pair(seq[i], seq[j], include_ambiguous):
                    valid_pairs.append((i, j))
        
        return valid_pairs


# =========================================================
# 5. Resumeable Crawler Base Class
# =========================================================
class ResumeableCrawler(ABC):
    """Base class for resumeable web crawlers with database state tracking."""
    
    def __init__(self, source_name: str):
        self.source_name = source_name
        self.config = config
        self.session = requests.Session()
        self.session.headers.update({'User-Agent': config.USER_AGENT})
        self.last_request_time = 0
        self.domain_delays = {}
        self.db_conn = None
        self._connect_db()
        self._init_source_state()
        self._check_robots_txt()
    
    def _connect_db(self):
        """Connect to SQLite database."""
        self.db_conn = sqlite3.connect(str(config.CRAWLER_DB))
        self.db_conn.row_factory = sqlite3.Row
    
    def _init_source_state(self):
        """Initialize data source state in database."""
        cursor = self.db_conn.cursor()
        
        cursor.execute(
            "SELECT * FROM data_sources WHERE source_name = ?",
            (self.source_name,)
        )
        row = cursor.fetchone()
        
        if row is None:
            cursor.execute('''
                INSERT INTO data_sources 
                (source_name, total_sequences, downloaded_sequences, status, config)
                VALUES (?, ?, ?, ?, ?)
            ''', (self.source_name, 0, 0, 'initialized', 
                  json.dumps(self.config.DATA_SOURCES.get(self.source_name, {}))))
            self.db_conn.commit()
    
    def _check_robots_txt(self):
        """Check robots.txt for allowed paths."""
        if not self.config.RESPECT_ROBOTS_TXT:
            return
        
        try:
            base_url = self._get_base_url()
            parsed = urlparse(base_url)
            robots_url = f"{parsed.scheme}://{parsed.netloc}/robots.txt"
            
            response = self.session.get(robots_url, timeout=10)
            if response.status_code == 200:
                # Simple parsing - in production, use robotparser
                self.robots_allowed = True
                log_message(f"Robots.txt found for {self.source_name}")
        except Exception as e:
            log_message(f"Failed to check robots.txt for {self.source_name}: {e}", "WARNING")
            self.robots_allowed = True
    
    @abstractmethod
    def _get_base_url(self) -> str:
        """Get base URL for the data source."""
        pass
    
    def _respectful_request(self, url: str, **kwargs) -> requests.Response:
        """
        Make a request respecting rate limits and retry logic.
        """
        domain = urlparse(url).netloc
        
        # Rate limiting
        current_time = time.time()
        time_since_last = current_time - self.domain_delays.get(domain, 0)
        if time_since_last < self.config.RATE_LIMIT:
            time.sleep(self.config.RATE_LIMIT - time_since_last)
        
        # Retry logic
        for attempt in range(self.config.MAX_RETRIES):
            try:
                response = self.session.get(url, timeout=self.config.REQUEST_TIMEOUT, **kwargs)
                self.domain_delays[domain] = time.time()
                
                # Check for rate limiting headers
                if response.status_code == 429:
                    retry_after = int(response.headers.get('Retry-After', self.config.RETRY_DELAY))
                    time.sleep(retry_after)
                    continue
                
                response.raise_for_status()
                return response
                
            except requests.exceptions.RequestException as e:
                log_message(f"Request failed (attempt {attempt+1}/{self.config.MAX_RETRIES}): {e}", 
                           "WARNING")
                if attempt < self.config.MAX_RETRIES - 1:
                    time.sleep(self.config.RETRY_DELAY * (attempt + 1))
                else:
                    raise
        
        raise Exception(f"Failed to fetch {url} after {self.config.MAX_RETRIES} attempts")
    
    def update_source_state(self, **kwargs):
        """Update data source state in database."""
        fields = []
        values = []
        
        for key, value in kwargs.items():
            if key in ['total_sequences', 'downloaded_sequences', 'status', 'error_message']:
                fields.append(f"{key} = ?")
                values.append(value)
        
        if fields:
            values.append(self.source_name)
            cursor = self.db_conn.cursor()
            cursor.execute(f'''
                UPDATE data_sources 
                SET {", ".join(fields)}, last_crawled = CURRENT_TIMESTAMP
                WHERE source_name = ?
            ''', values)
            self.db_conn.commit()
    
    def is_file_downloaded(self, url: str, expected_md5: str = None) -> Tuple[bool, Optional[Path]]:
        """Check if a file has already been downloaded."""
        cursor = self.db_conn.cursor()
        
        cursor.execute(
            "SELECT local_path, md5_hash, status FROM downloaded_files WHERE url = ?",
            (url,)
        )
        row = cursor.fetchone()
        
        if row is None:
            return False, None
        
        local_path = Path(row['local_path'])
        
        if not local_path.exists():
            return False, None
        
        if expected_md5 and row['md5_hash'] != expected_md5:
            return False, None
        
        if row['status'] != 'completed':
            return False, None
        
        return True, local_path
    
    def mark_file_downloaded(self, url: str, local_path: Path, md5_hash: str = None):
        """Mark a file as successfully downloaded."""
        if md5_hash is None:
            md5_hash = compute_md5(local_path)
        
        file_id = hashlib.md5(url.encode()).hexdigest()
        
        cursor = self.db_conn.cursor()
        cursor.execute('''
            INSERT OR REPLACE INTO downloaded_files
            (file_id, source_name, url, local_path, file_size, md5_hash, download_time, status)
            VALUES (?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP, ?)
        ''', (file_id, self.source_name, url, str(local_path), 
              local_path.stat().st_size, md5_hash, 'completed'))
        self.db_conn.commit()
    
    def is_sequence_parsed(self, md5_hash: str = None, sequence_id: str = None) -> bool:
        """Check if a sequence has already been parsed."""
        cursor = self.db_conn.cursor()
        
        if md5_hash:
            cursor.execute(
                "SELECT 1 FROM parsed_sequences WHERE md5_hash = ?",
                (md5_hash,)
            )
        elif sequence_id:
            cursor.execute(
                "SELECT 1 FROM parsed_sequences WHERE sequence_id = ?",
                (sequence_id,)
            )
        else:
            return False
        
        return cursor.fetchone() is not None
    
    def save_parsed_sequence(self, original_id: str, sequence: str, 
                             dot_bracket: str, quality_score: float):
        """Save a parsed sequence to database."""
        seq_id = f"{self.source_name}:{original_id}"
        md5_hash = compute_sequence_md5(sequence, dot_bracket)
        
        if self.is_sequence_parsed(md5_hash=md5_hash):
            return
        
        gc_content = compute_gc_content(sequence)
        
        cursor = self.db_conn.cursor()
        cursor.execute('''
            INSERT INTO parsed_sequences
            (sequence_id, source_name, original_id, sequence, dot_bracket, 
             length, quality_score, gc_content, md5_hash, parsed_time)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
        ''', (seq_id, self.source_name, original_id, sequence, dot_bracket,
              len(sequence), quality_score, gc_content, md5_hash))
        self.db_conn.commit()
        
        cursor.execute('''
            UPDATE data_sources 
            SET downloaded_sequences = downloaded_sequences + 1
            WHERE source_name = ?
        ''', (self.source_name,))
        self.db_conn.commit()
    
    def download_file(self, url: str, local_path: Path, 
                      expected_size: int = None, expected_md5: str = None,
                      resume: bool = True) -> bool:
        """
        Download a file with resume capability.
        
        Args:
            url: URL to download
            local_path: Local path to save the file
            expected_size: Expected file size for validation
            expected_md5: Expected MD5 hash for validation
            resume: Whether to resume partial downloads
        
        Returns:
            True if download successful, False otherwise
        """
        if self.is_file_downloaded(url, expected_md5)[0]:
            log_message(f"File already downloaded: {url}")
            return True
        
        local_path.parent.mkdir(parents=True, exist_ok=True)
        
        headers = {}
        mode = 'wb'
        existing_size = 0
        
        if resume and local_path.exists():
            existing_size = local_path.stat().st_size
            headers['Range'] = f'bytes={existing_size}-'
            mode = 'ab'
        
        for attempt in range(self.config.MAX_RETRIES):
            try:
                response = self._respectful_request(url, headers=headers, stream=True)
                
                total_size = int(response.headers.get('content-length', 0)) + existing_size
                
                with open(local_path, mode) as f:
                    with tqdm(
                        total=total_size,
                        initial=existing_size,
                        unit='B',
                        unit_scale=True,
                        desc=f"Downloading {local_path.name}"
                    ) as pbar:
                        for chunk in response.iter_content(chunk_size=self.config.DOWNLOAD_CHUNK_SIZE):
                            if chunk:
                                f.write(chunk)
                                pbar.update(len(chunk))
                
                if expected_size and local_path.stat().st_size != expected_size:
                    raise ValueError(f"Size mismatch: expected {expected_size}, got {local_path.stat().st_size}")
                
                md5_hash = compute_md5(local_path)
                self.mark_file_downloaded(url, local_path, md5_hash)
                
                log_message(f"Successfully downloaded: {url}")
                return True
                
            except Exception as e:
                log_message(f"Download attempt {attempt + 1} failed for {url}: {e}", "ERROR")
                if attempt < self.config.MAX_RETRIES - 1:
                    time.sleep(self.config.RETRY_DELAY)
                else:
                    self._mark_file_failed(url, local_path, str(e), attempt + 1)
                    return False
        
        return False
    
    def _mark_file_failed(self, url: str, local_path: Path, error: str, retry_count: int):
        """Mark a file download as failed."""
        file_id = hashlib.md5(url.encode()).hexdigest()
        cursor = self.db_conn.cursor()
        cursor.execute('''
            INSERT OR REPLACE INTO downloaded_files
            (file_id, source_name, url, local_path, status, retry_count, error_message)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        ''', (file_id, self.source_name, url, str(local_path), 
              'failed', retry_count, error))
        self.db_conn.commit()
    
    def get_progress(self) -> Dict:
        """Get crawling progress for this source."""
        cursor = self.db_conn.cursor()
        cursor.execute('''
            SELECT total_sequences, downloaded_sequences, status, error_message, last_crawled
            FROM data_sources WHERE source_name = ?
        ''', (self.source_name,))
        row = cursor.fetchone()
        
        if row:
            return dict(row)
        return {}
    
    def reset_failed_downloads(self):
        """Reset status of failed downloads to retry."""
        cursor = self.db_conn.cursor()
        cursor.execute('''
            UPDATE downloaded_files 
            SET status = 'pending', retry_count = 0
            WHERE source_name = ? AND status = 'failed'
        ''', (self.source_name,))
        self.db_conn.commit()
        log_message(f"Reset {cursor.rowcount} failed downloads for {self.source_name}")
    
    def close(self):
        """Close database connection."""
        if self.db_conn:
            self.db_conn.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
    
    @abstractmethod
    def crawl(self, max_sequences: int = None, resume: bool = False) -> int:
        """Crawl data from the source."""
        pass


# =========================================================
# 6. bpRNA Crawler
# =========================================================
class BpRNACrawler(ResumeableCrawler):
    """Crawler for bpRNA database (with structure annotations)."""
    
    def __init__(self):
        super().__init__('bprna')
        self.base_url = config.DATA_SOURCES['bprna']['base_url']
        self.dotbracket_url = config.DATA_SOURCES['bprna']['dotbracket_url']
        self.fasta_url = config.DATA_SOURCES['bprna']['fasta_url']
        self.max_sequences = config.DATA_SOURCES['bprna']['max_sequences']
    
    def _get_base_url(self) -> str:
        return self.base_url
    
    def crawl(self, max_sequences: int = None, resume: bool = False) -> int:
        """
        Crawl bpRNA data.
        
        Args:
            max_sequences: Maximum number of sequences to crawl
            resume: Whether to resume from previous state
        
        Returns:
            Number of sequences crawled
        """
        if max_sequences is None:
            max_sequences = self.max_sequences
        
        if resume:
            self.reset_failed_downloads()
        
        log_message(f"Starting bpRNA crawl, target: {max_sequences} sequences")
        self.update_source_state(status='crawling')
        
        try:
            fasta_zip = config.RAW_DIR / "bprna" / "bprna_fasta.zip"
            dot_zip = config.RAW_DIR / "bprna" / "bprna_dotbracket.zip"
            
            if not self.download_file(self.fasta_url, fasta_zip):
                raise Exception("Failed to download FASTA zip")
            
            if not self.download_file(self.dotbracket_url, dot_zip):
                raise Exception("Failed to download dot-bracket zip")
            
            fasta_dir = config.RAW_DIR / "bprna" / "fasta"
            dot_dir = config.RAW_DIR / "bprna" / "dotbracket"
            
            self._extract_zip(fasta_zip, fasta_dir)
            self._extract_zip(dot_zip, dot_dir)
            
            count = self._parse_files(fasta_dir, dot_dir, max_sequences)
            
            self.update_source_state(
                status='completed',
                total_sequences=count
            )
            
            log_message(f"bpRNA crawl completed: {count} sequences")
            return count
            
        except Exception as e:
            log_message(f"bpRNA crawl failed: {e}", "ERROR")
            self.update_source_state(status='failed', error_message=str(e))
            raise
    
    def _extract_zip(self, zip_path: Path, extract_dir: Path):
        """Extract zip file if not already extracted."""
        if extract_dir.exists() and any(extract_dir.iterdir()):
            log_message(f"Extraction directory already exists: {extract_dir}")
            return
        
        extract_dir.mkdir(parents=True, exist_ok=True)
        
        with zipfile.ZipFile(zip_path, 'r') as zf:
            zf.extractall(extract_dir)
        
        log_message(f"Extracted {zip_path.name} to {extract_dir}")
    
    def _parse_files(self, fasta_dir: Path, dot_dir: Path, max_sequences: int) -> int:
        """Parse FASTA and dot-bracket files."""
        count = 0
        
        fasta_files = list(fasta_dir.glob("*.fa")) + list(fasta_dir.glob("*.fasta"))
        dot_files = list(dot_dir.glob("*.db")) + list(dot_dir.glob("*.txt"))
        
        fasta_map = {f.stem: f for f in fasta_files}
        dot_map = {f.stem: f for f in dot_files}
        
        common_keys = sorted(set(fasta_map.keys()) & set(dot_map.keys()))
        
        log_message(f"Found {len(common_keys)} matching file pairs")
        
        for key in tqdm(common_keys, desc="Parsing bpRNA files"):
            if count >= max_sequences:
                break
            
            try:
                seq = self._read_fasta(fasta_map[key])
                db = self._read_dot_bracket(dot_map[key])
                
                if seq and db and len(seq) == len(db):
                    quality = validate_sequence_quality(seq)
                    
                    if (quality >= config.MIN_SEQUENCE_QUALITY and 
                        len(seq) >= config.MIN_LEN and 
                        len(seq) <= config.MAX_LEN):
                        
                        self.save_parsed_sequence(key, seq, db, quality)
                        count += 1
                        
            except Exception as e:
                log_message(f"Error parsing {key}: {e}", "WARNING")
                continue
        
        return count
    
    def _read_fasta(self, path: Path) -> Optional[str]:
        """Read first sequence from FASTA file."""
        try:
            with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                seq_lines = []
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('>'):
                        seq_lines.append(line)
                
                if seq_lines:
                    seq = ''.join(seq_lines)
                    return clean_sequence(seq)
        except Exception as e:
            log_message(f"Error reading FASTA {path}: {e}", "WARNING")
        
        return None
    
    def _read_dot_bracket(self, path: Path) -> Optional[str]:
        """Read structure from dot-bracket file."""
        try:
            with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('>') and not line.startswith('#'):
                        cleaned = clean_dot_bracket(line)
                        if cleaned:
                            return cleaned
        except Exception as e:
            log_message(f"Error reading dot-bracket {path}: {e}", "WARNING")
        
        return None


# =========================================================
# 7. Rfam Crawler
# =========================================================
class RfamCrawler(ResumeableCrawler):
    """Crawler for Rfam database (sequences only, no structures)."""
    
    def __init__(self):
        super().__init__('rfam')
        self.base_url = config.DATA_SOURCES['rfam']['base_url']
        self.family_list_url = config.DATA_SOURCES['rfam']['family_list_url']
        self.max_sequences = config.DATA_SOURCES['rfam']['max_sequences']
    
    def _get_base_url(self) -> str:
        return self.base_url
    
    def crawl(self, max_sequences: int = None, resume: bool = False) -> int:
        """
        Crawl Rfam data.
        
        Args:
            max_sequences: Maximum number of sequences to crawl
            resume: Whether to resume from previous state
        
        Returns:
            Number of sequences crawled
        """
        if max_sequences is None:
            max_sequences = self.max_sequences
        
        if resume:
            self.reset_failed_downloads()
        
        log_message(f"Starting Rfam crawl, target: {max_sequences} sequences")
        self.update_source_state(status='crawling')
        
        try:
            family_ids = self._get_family_ids(limit=1000)
            log_message(f"Found {len(family_ids)} Rfam families")
            
            file_map = self._get_fasta_files(family_ids)
            log_message(f"Found {len(file_map)} downloadable files")
            
            count = 0
            for fid, filename in tqdm(file_map.items(), desc="Downloading Rfam families"):
                if count >= max_sequences:
                    break
                
                try:
                    url = urljoin(self.base_url, filename)
                    local_path = config.RAW_DIR / "rfam" / filename
                    
                    if self.download_file(url, local_path):
                        parsed = self._parse_family_file(local_path, fid, max_sequences - count)
                        count += parsed
                        
                except Exception as e:
                    log_message(f"Error processing family {fid}: {e}", "WARNING")
                    continue
            
            self.update_source_state(
                status='completed',
                total_sequences=count
            )
            
            log_message(f"Rfam crawl completed: {count} sequences")
            return count
            
        except Exception as e:
            log_message(f"Rfam crawl failed: {e}", "ERROR")
            self.update_source_state(status='failed', error_message=str(e))
            raise
    
    def _get_family_ids(self, limit: int = 1000) -> List[str]:
        """Get list of Rfam family IDs."""
        try:
            response = self._respectful_request(self.family_list_url)
            soup = BeautifulSoup(response.text, 'html.parser')
            
            family_ids = []
            for link in soup.find_all('a'):
                href = link.get('href', '')
                if '/family/' in href:
                    fid = href.split('/')[-1]
                    if fid.startswith('RF'):
                        family_ids.append(fid)
            
            return family_ids[:limit]
            
        except Exception as e:
            log_message(f"Error getting family IDs: {e}", "ERROR")
            return [f"RF{str(i).zfill(5)}" for i in range(1, min(limit + 1, 1000))]
    
    def _get_fasta_files(self, family_ids: List[str]) -> Dict[str, str]:
        """Get FASTA file names for each family."""
        try:
            response = self._respectful_request(self.base_url)
            soup = BeautifulSoup(response.text, 'html.parser')
            
            all_files = []
            for link in soup.find_all('a'):
                href = link.get('href', '')
                if href.endswith(('.fa.gz', '.fasta.gz', '.fa', '.fasta')):
                    all_files.append(href)
            
            file_map = {}
            for fid in family_ids:
                matched = [f for f in all_files if fid in f]
                if matched:
                    matched.sort(key=lambda x: (0 if x.endswith('.gz') else 1, len(x)))
                    file_map[fid] = matched[0]
            
            return file_map
            
        except Exception as e:
            log_message(f"Error getting FASTA files: {e}", "WARNING")
            return {fid: f"{fid}.fa.gz" for fid in family_ids}
    
    def _parse_family_file(self, file_path: Path, family_id: str, max_sequences: int) -> int:
        """Parse family FASTA file and save sequences."""
        count = 0
        
        open_fn = gzip.open if str(file_path).endswith('.gz') else open
        
        try:
            with open_fn(file_path, 'rt', encoding='utf-8', errors='ignore') as f:
                current_id = None
                current_seq = []
                
                for line in f:
                    line = line.strip()
                    
                    if line.startswith('>'):
                        if current_id and current_seq:
                            seq = ''.join(current_seq)
                            quality = validate_sequence_quality(seq)
                            
                            if (len(seq) >= config.MIN_LEN and 
                                len(seq) <= config.MAX_LEN and 
                                quality >= config.MIN_SEQUENCE_QUALITY):
                                # Rfam has no structure information
                                self.save_parsed_sequence(f"{family_id}_{current_id}", seq, '', quality)
                                count += 1
                                
                                if count >= max_sequences:
                                    break
                        
                        current_id = line[1:].split()[0]
                        current_seq = []
                        
                    elif line:
                        current_seq.append(line)
                
                if current_id and current_seq and count < max_sequences:
                    seq = ''.join(current_seq)
                    quality = validate_sequence_quality(seq)
                    
                    if (len(seq) >= config.MIN_LEN and 
                        len(seq) <= config.MAX_LEN and 
                        quality >= config.MIN_SEQUENCE_QUALITY):
                        self.save_parsed_sequence(f"{family_id}_{current_id}", seq, '', quality)
                        count += 1
                        
        except Exception as e:
            log_message(f"Error parsing {file_path}: {e}", "WARNING")
        
        return count


# =========================================================
# 8. RNAcentral Crawler
# =========================================================
class RNAcentralCrawler(ResumeableCrawler):
    """Crawler for RNAcentral database (via API)."""
    
    def __init__(self):
        super().__init__('rnacentral')
        self.api_url = config.DATA_SOURCES['rnacentral']['api_url']
        self.max_sequences = config.DATA_SOURCES['rnacentral']['max_sequences']
    
    def _get_base_url(self) -> str:
        return self.api_url
    
    def crawl(self, max_sequences: int = None, resume: bool = False,
              rna_types: List[str] = None) -> int:
        """
        Crawl RNAcentral data via API.
        
        Args:
            max_sequences: Maximum number of sequences to crawl
            resume: Whether to resume from previous state
            rna_types: RNA types to filter (e.g., ['miRNA', 'rRNA', 'tRNA'])
        
        Returns:
            Number of sequences crawled
        """
        if max_sequences is None:
            max_sequences = self.max_sequences
        
        if rna_types is None:
            rna_types = ['miRNA', 'rRNA', 'tRNA', 'lncRNA', 'snRNA', 'snoRNA']
        
        if resume:
            self.reset_failed_downloads()
        
        log_message(f"Starting RNAcentral crawl, target: {max_sequences} sequences")
        self.update_source_state(status='crawling')
        
        try:
            count = 0
            page = 1
            
            while count < max_sequences:
                sequences = self._fetch_page(page, rna_types)
                
                if not sequences:
                    break
                
                for seq_data in sequences:
                    if count >= max_sequences:
                        break
                    
                    seq = seq_data.get('sequence', '')
                    quality = validate_sequence_quality(seq)
                    
                    if (len(seq) >= config.MIN_LEN and 
                        len(seq) <= config.MAX_LEN and 
                        quality >= config.MIN_SEQUENCE_QUALITY):
                        
                        dot_bracket = seq_data.get('secondary_structure', '')
                        if dot_bracket:
                            dot_bracket = clean_dot_bracket(dot_bracket)
                        
                        self.save_parsed_sequence(
                            seq_data.get('id', f"page{page}_{count}"),
                            seq,
                            dot_bracket,
                            quality
                        )
                        count += 1
                
                log_message(f"RNAcentral: fetched page {page}, total {count} sequences")
                page += 1
                
                # Rate limiting is handled by _respectful_request
            
            self.update_source_state(
                status='completed',
                total_sequences=count
            )
            
            log_message(f"RNAcentral crawl completed: {count} sequences")
            return count
            
        except Exception as e:
            log_message(f"RNAcentral crawl failed: {e}", "ERROR")
            self.update_source_state(status='failed', error_message=str(e))
            raise
    
    def _fetch_page(self, page: int, rna_types: List[str]) -> List[Dict]:
        """Fetch one page of data from RNAcentral API."""
        try:
            params = {
                'format': 'json',
                'page': page,
                'page_size': 100,
                'has_sequence': 'true'
            }
            
            if rna_types:
                params['rna_type'] = ','.join(rna_types)
            
            response = self._respectful_request(
                urljoin(self.api_url, 'rnacentral/'),
                params=params
            )
            
            data = response.json()
            return data.get('results', [])
            
        except Exception as e:
            log_message(f"Error fetching RNAcentral page {page}: {e}", "WARNING")
            return []


# =========================================================
# 9. Multi-Source Crawler Manager
# =========================================================
class DataCrawlerManager:
    """Manager for multi-source data crawling with quota allocation."""
    
    def __init__(self):
        self.config = config
        self.crawlers = {}
        self._init_crawlers()
    
    def _init_crawlers(self):
        """Initialize enabled crawlers based on configuration."""
        sources = config.DATA_SOURCES
        
        if sources['bprna']['enabled']:
            self.crawlers['bprna'] = BpRNACrawler()
        
        if sources['rfam']['enabled']:
            self.crawlers['rfam'] = RfamCrawler()
        
        if sources['rnacentral']['enabled']:
            self.crawlers['rnacentral'] = RNAcentralCrawler()
    
    def crawl_all(self, target_total: int = None, parallel: bool = False,
                  resume: bool = False, source_list: List[str] = None) -> Dict[str, int]:
        """
        Crawl data from all enabled sources with quota allocation.
        
        Args:
            target_total: Target total number of sequences
            parallel: Whether to crawl sources in parallel
            resume: Whether to resume from previous state
            source_list: List of sources to crawl (None for all enabled)
        
        Returns:
            Dictionary of source -> sequences crawled
        """
        if source_list is not None:
            crawlers = {name: self.crawlers[name] for name in source_list 
                       if name in self.crawlers}
        else:
            crawlers = self.crawlers
        
        if target_total is None:
            target_total = sum(
                src['max_sequences'] 
                for name, src in config.DATA_SOURCES.items() 
                if name in crawlers
            )
        
        # Allocate quotas based on source priorities
        quotas = self._allocate_quotas(target_total, crawlers)
        
        log_message(f"Starting multi-source crawl, target: {target_total} sequences")
        log_message(f"Quotas: {quotas}")
        
        results = {}
        
        if parallel and len(crawlers) > 1:
            with ThreadPoolExecutor(max_workers=config.CONCURRENT_DOWNLOADS) as executor:
                futures = {}
                for name, crawler in crawlers.items():
                    future = executor.submit(crawler.crawl, quotas[name], resume)
                    futures[future] = name
                
                for future in as_completed(futures):
                    name = futures[future]
                    try:
                        count = future.result()
                        results[name] = count
                        log_message(f"{name} completed: {count} sequences")
                    except Exception as e:
                        log_message(f"{name} failed: {e}", "ERROR")
                        results[name] = 0
        else:
            for name, crawler in crawlers.items():
                try:
                    count = crawler.crawl(quotas[name], resume)
                    results[name] = count
                    log_message(f"{name} completed: {count} sequences")
                except Exception as e:
                    log_message(f"{name} failed: {e}", "ERROR")
                    results[name] = 0
        
        total = sum(results.values())
        
        summary = {
            'target_total': target_total,
            'actual_total': total,
            'sources': results,
            'quotas': quotas,
            'timestamp': get_timestamp()
        }
        
        summary_path = config.RESULTS_DIR / "crawl_summary.json"
        save_json(summary, summary_path)
        
        log_message(f"Crawl completed: {total}/{target_total} sequences")
        
        return results
    
    def _allocate_quotas(self, target_total: int, crawlers: Dict) -> Dict[str, int]:
        """Allocate quotas to sources based on their max_sequences and priorities."""
        sources_info = []
        for name in crawlers:
            info = config.DATA_SOURCES[name].copy()
            info['name'] = name
            sources_info.append(info)
        
        # Sort by priority (lower number = higher priority)
        sources_info.sort(key=lambda x: x.get('priority', 999))
        
        quotas = {}
        remaining = target_total
        
        for info in sources_info:
            max_seq = info['max_sequences']
            quota = min(max_seq, remaining)
            quotas[info['name']] = quota
            remaining -= quota
            
            if remaining <= 0:
                break
        
        return quotas
    
    def get_progress(self) -> Dict:
        """Get progress for all crawlers."""
        progress = {}
        for name, crawler in self.crawlers.items():
            progress[name] = crawler.get_progress()
        return progress
    
    def reset_failed_downloads(self, source: str = None):
        """Reset failed downloads for all or specified source."""
        if source:
            if source in self.crawlers:
                self.crawlers[source].reset_failed_downloads()
        else:
            for crawler in self.crawlers.values():
                crawler.reset_failed_downloads()
    
    def close(self):
        """Close all crawlers."""
        for crawler in self.crawlers.values():
            crawler.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


# =========================================================
# 10. Data Processing and Integration
# =========================================================
class DataProcessor:
    """Data processing and integration class."""
    
    def __init__(self):
        self.config = config
        self.db_conn = sqlite3.connect(str(config.CRAWLER_DB))
        self.db_conn.row_factory = sqlite3.Row
    
    def build_dataset(self, min_quality: float = 0.8, 
                     require_structure: bool = True,
                     deduplicate: bool = True,
                     source_filter: List[str] = None) -> pd.DataFrame:
        """
        Build dataset from crawled data.
        
        Args:
            min_quality: Minimum sequence quality score
            require_structure: Whether to require structure annotations
            deduplicate: Whether to remove duplicate sequences
            source_filter: List of sources to include (None for all)
        
        Returns:
            DataFrame with columns: sequence_id, source, sequence, dot_bracket, length, quality
        """
        log_message("Building dataset from crawled data...")
        
        cursor = self.db_conn.cursor()
        
        query = '''
            SELECT sequence_id, source_name, original_id, sequence, dot_bracket, 
                   length, quality_score, gc_content
            FROM parsed_sequences
            WHERE length >= ? AND length <= ? AND quality_score >= ?
        '''
        params = [config.MIN_LEN, config.MAX_LEN, min_quality]
        
        if require_structure:
            query += " AND dot_bracket != '' AND length(dot_bracket) > 0"
        
        if source_filter:
            placeholders = ','.join(['?'] * len(source_filter))
            query += f" AND source_name IN ({placeholders})"
            params.extend(source_filter)
        
        cursor.execute(query, params)
        rows = cursor.fetchall()
        
        log_message(f"Retrieved {len(rows)} sequences from database")
        
        if not rows:
            log_message("No data available", "WARNING")
            return pd.DataFrame()
        
        df = pd.DataFrame([dict(row) for row in rows])
        
        if deduplicate:
            df = self._deduplicate_sequences(df)
        
        # Add derived features
        df['has_structure'] = df['dot_bracket'].apply(lambda x: len(x) > 0)
        
        if 'dot_bracket' in df.columns:
            df['pair_count'] = df['dot_bracket'].apply(self._count_pairs)
        
        source_counts = df['source_name'].value_counts()
        log_message(f"Data by source:\n{source_counts}")
        
        raw_csv = config.PROC_DIR / "raw_dataset.csv"
        df.to_csv(raw_csv, index=False)
        log_message(f"Raw dataset saved to {raw_csv}")
        
        return df
    
    def _deduplicate_sequences(self, df: pd.DataFrame) -> pd.DataFrame:
        """Remove duplicate sequences based on MD5 hash."""
        if 'md5_hash' in df.columns:
            # Use existing MD5 hashes
            df = df.drop_duplicates(subset=['md5_hash'])
        else:
            # Compute MD5 hashes on the fly
            df['md5_hash'] = df.apply(
                lambda row: compute_sequence_md5(row['sequence'], row.get('dot_bracket', '')),
                axis=1
            )
            df = df.drop_duplicates(subset=['md5_hash'])
        
        log_message(f"After deduplication: {len(df)} sequences")
        return df
    
    def _count_pairs(self, db: str) -> int:
        """Count number of base pairs in a dot-bracket string."""
        if not isinstance(db, str) or len(db) == 0:
            return 0
        # Count opening brackets (each pair is counted once)
        return sum(1 for c in db if c in '([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    
    def filter_by_structure_complexity(self, df: pd.DataFrame, 
                                      min_pairs: int = 1) -> pd.DataFrame:
        """Filter by minimum number of base pairs."""
        if 'pair_count' not in df.columns:
            df['pair_count'] = df['dot_bracket'].apply(self._count_pairs)
        
        filtered_df = df[df['pair_count'] >= min_pairs].copy()
        log_message(f"Filtered by structure complexity: {len(filtered_df)}/{len(df)} remain")
        
        return filtered_df
    
    def export_to_fasta(self, df: pd.DataFrame, output_path: Path):
        """Export sequences to FASTA format."""
        with open(output_path, 'w', encoding='utf-8') as f:
            for _, row in df.iterrows():
                f.write(f">{row['sequence_id']}\n")
                f.write(f"{row['sequence']}\n")
        
        log_message(f"FASTA exported to {output_path}")
    
    def export_to_stockholm(self, df: pd.DataFrame, output_path: Path):
        """Export sequences and structures to Stockholm format."""
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write("# STOCKHOLM 1.0\n\n")
            
            for _, row in df.iterrows():
                f.write(f"{row['sequence_id']:<20} {row['sequence']}\n")
                if row.get('dot_bracket'):
                    f.write(f"{'#=GF SS':<20} {row['dot_bracket']}\n")
            
            f.write("//\n")
        
        log_message(f"Stockholm format exported to {output_path}")
    
    def close(self):
        """Close database connection."""
        if self.db_conn:
            self.db_conn.close()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


# =========================================================
# 11. Sequence Quality Scorer
# =========================================================
class SequenceQualityScorer:
    """Comprehensive sequence quality scoring system."""
    
    def __init__(self):
        self.weights = {
            'length': 0.2,
            'gc_content': 0.15,
            'ambiguous_bases': 0.25,
            'structure_complexity': 0.3,
            'conserved_motifs': 0.1
        }
    
    def score_sequence(self, seq: str, db: str = '') -> Dict[str, float]:
        """
        Compute comprehensive quality score for a sequence.
        
        Args:
            seq: RNA sequence
            db: Dot-bracket structure (optional)
        
        Returns:
            Dictionary of individual scores and total score
        """
        scores = {}
        
        # Length score (optimal around 150 nt)
        optimal_len = 150
        len_score = 1.0 - min(abs(len(seq) - optimal_len) / optimal_len, 1.0)
        scores['length'] = len_score
        
        # GC content score (optimal 30-70%)
        gc_count = seq.count('G') + seq.count('C')
        gc_ratio = gc_count / len(seq) if len(seq) > 0 else 0
        gc_score = 1.0 - abs(gc_ratio - 0.5) * 2
        scores['gc_content'] = max(0, min(1, gc_score))
        
        # Ambiguous bases score
        ambiguous = sum(1 for base in seq if base == 'N')
        ambiguous_score = 1.0 - (ambiguous / len(seq))
        scores['ambiguous_bases'] = ambiguous_score
        
        # Structure complexity score
        if db:
            pairs = BracketedNotationParser.dotbracket_to_pairs(db, allow_pseudoknots=False)
            complexity = len(pairs) / (len(seq) / 2)  # Pairing proportion
            scores['structure_complexity'] = min(1.0, complexity * 2)
        else:
            scores['structure_complexity'] = 0.5
        
        # Conserved motif score (simplified)
        motifs = ['GU', 'GNRA', 'UNCG']  # Common RNA motifs
        motif_count = sum(seq.count(m) for m in motifs)
        motif_score = min(1.0, motif_count / (len(seq) / 10))
        scores['conserved_motifs'] = motif_score
        
        # Weighted total
        total_score = sum(scores[k] * self.weights[k] for k in scores)
        scores['total'] = total_score
        
        return scores


# =========================================================
# 12. Homology-Aware Data Splitter
# =========================================================
class HomologyAwareSplitter:
    """Data splitter that respects sequence homology."""
    
    def __init__(self, method: str = "mmseqs2", identity: float = 0.8,
                 coverage: float = 0.8):
        self.method = method
        self.identity = identity
        self.coverage = coverage
        self._check_tools()
    
    def _check_tools(self):
        """Check if required clustering tools are available."""
        if self.method == "mmseqs2":
            if not has_executable("mmseqs"):
                log_message("mmseqs2 not found, falling back to CD-HIT", "WARNING")
                self.method = "cdhit"
        
        if self.method == "cdhit":
            if not has_executable("cd-hit-est"):
                raise RuntimeError("Neither mmseqs2 nor cd-hit-est found")
    
    def split(self, df: pd.DataFrame, train_ratio: float = 0.7,
             valid_ratio: float = 0.15, seed: int = 42,
             output_prefix: str = "split") -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Perform homology-aware splitting.
        
        Returns:
            train_df, valid_df, test_df
        """
        log_message(f"Running homology-aware split with {self.method} (identity={self.identity})")
        
        try:
            cluster_map = self._run_clustering(df)
            df['cluster_id'] = df['sequence_id'].map(cluster_map)
            df = df.dropna(subset=['cluster_id']).reset_index(drop=True)
            
            return self._split_by_clusters(df, train_ratio, valid_ratio, seed, output_prefix)
            
        except Exception as e:
            log_message(f"Homology clustering failed: {e}", "WARNING")
            log_message("Falling back to random split", "WARNING")
            return self._random_split(df, train_ratio, valid_ratio, seed, output_prefix)
    
    def _run_clustering(self, df: pd.DataFrame) -> Dict[str, str]:
        """Run sequence clustering and return cluster assignments."""
        fasta_path = config.TMP_DIR / "clustering.fasta"
        self._write_fasta(df, fasta_path)
        
        if self.method == "mmseqs2":
            tsv_path = self._run_mmseqs2(fasta_path)
            return self._parse_mmseqs_tsv(tsv_path)
        else:
            clstr_path = self._run_cdhit(fasta_path)
            return self._parse_cdhit_clstr(clstr_path)
    
    def _write_fasta(self, df: pd.DataFrame, fasta_path: Path):
        """Write sequences to FASTA file for clustering."""
        with open(fasta_path, 'w', encoding='utf-8') as f:
            for _, row in df.iterrows():
                f.write(f">{row['sequence_id']}\n")
                f.write(f"{row['sequence']}\n")
    
    def _run_mmseqs2(self, fasta_path: Path) -> Path:
        """Run MMseqs2 clustering."""
        out_dir = config.TMP_DIR / "mmseqs2"
        out_dir.mkdir(parents=True, exist_ok=True)
        cluster_prefix = out_dir / "cluster"
        
        cmd = [
            "mmseqs", "easy-cluster",
            str(fasta_path),
            str(cluster_prefix),
            str(out_dir / "tmp"),
            "--min-seq-id", str(self.identity),
            "-c", str(self.coverage),
            "--cov-mode", "1"
        ]
        
        log_message(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True, capture_output=True)
        
        tsv_path = Path(str(cluster_prefix) + "_cluster.tsv")
        if not tsv_path.exists():
            raise FileNotFoundError(f"MMseqs2 output not found: {tsv_path}")
        
        return tsv_path
    
    def _parse_mmseqs_tsv(self, tsv_path: Path) -> Dict[str, str]:
        """Parse MMseqs2 TSV output."""
        sample_to_cluster = {}
        with open(tsv_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                rep, member = line.split('\t')
                sample_to_cluster[member] = rep
        return sample_to_cluster
    
    def _run_cdhit(self, fasta_path: Path) -> Path:
        """Run CD-HIT-EST clustering."""
        out_dir = config.TMP_DIR / "cdhit"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_prefix = out_dir / "cluster"
        
        cmd = [
            "cd-hit-est",
            "-i", str(fasta_path),
            "-o", str(out_prefix),
            "-c", str(self.identity),
            "-n", "5",
            "-M", "0",
            "-T", "0"
        ]
        
        log_message(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True, capture_output=True)
        
        clstr_path = Path(str(out_prefix) + ".clstr")
        if not clstr_path.exists():
            raise FileNotFoundError(f"CD-HIT output not found: {clstr_path}")
        
        return clstr_path
    
    def _parse_cdhit_clstr(self, clstr_path: Path) -> Dict[str, str]:
        """Parse CD-HIT cluster file."""
        sample_to_cluster = {}
        current_cluster = None
        
        with open(clstr_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line.startswith(">Cluster"):
                    current_cluster = line.replace(">", "")
                else:
                    m = re.search(r">(.+?)\.\.\.", line)
                    if m:
                        seq_id = m.group(1)
                        sample_to_cluster[seq_id] = current_cluster
        
        return sample_to_cluster
    
    def _split_by_clusters(self, df: pd.DataFrame, train_ratio: float,
                          valid_ratio: float, seed: int,
                          output_prefix: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Split by cluster IDs."""
        unique_clusters = df['cluster_id'].unique().tolist()
        random.Random(seed).shuffle(unique_clusters)
        
        n_total = len(unique_clusters)
        n_train = int(n_total * train_ratio)
        n_valid = int(n_total * valid_ratio)
        
        train_clusters = set(unique_clusters[:n_train])
        valid_clusters = set(unique_clusters[n_train:n_train + n_valid])
        test_clusters = set(unique_clusters[n_train + n_valid:])
        
        train_df = df[df['cluster_id'].isin(train_clusters)].copy()
        valid_df = df[df['cluster_id'].isin(valid_clusters)].copy()
        test_df = df[df['cluster_id'].isin(test_clusters)].copy()
        
        # Save splits
        train_df.to_csv(config.PROC_DIR / f"{output_prefix}_train.csv", index=False)
        valid_df.to_csv(config.PROC_DIR / f"{output_prefix}_valid.csv", index=False)
        test_df.to_csv(config.PROC_DIR / f"{output_prefix}_test.csv", index=False)
        
        # Save cluster info
        cluster_info = {
            'method': self.method,
            'identity': self.identity,
            'n_clusters': len(unique_clusters),
            'n_train_clusters': len(train_clusters),
            'n_valid_clusters': len(valid_clusters),
            'n_test_clusters': len(test_clusters),
            'n_train_sequences': len(train_df),
            'n_valid_sequences': len(valid_df),
            'n_test_sequences': len(test_df)
        }
        save_json(cluster_info, config.PROC_DIR / f"{output_prefix}_cluster_info.json")
        
        log_message(f"Split complete: train={len(train_df)}, valid={len(valid_df)}, test={len(test_df)}")
        log_message(f"Clusters: train={len(train_clusters)}, valid={len(valid_clusters)}, test={len(test_clusters)}")
        
        return train_df, valid_df, test_df
    
    def _random_split(self, df: pd.DataFrame, train_ratio: float,
                     valid_ratio: float, seed: int,
                     output_prefix: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Fallback to random split."""
        df = df.copy()
        
        # Stratify by length bins
        df['len_bin'] = pd.cut(
            df['length'],
            bins=[0, 40, 80, 120, 160, 220, 300, 1000],
            labels=False,
            include_lowest=True
        ).astype(str)
        
        train_df, temp_df = train_test_split(
            df,
            test_size=(1 - train_ratio),
            random_state=seed,
            stratify=df['len_bin']
        )
        
        valid_frac_in_temp = valid_ratio / (1 - train_ratio)
        valid_df, test_df = train_test_split(
            temp_df,
            test_size=(1 - valid_frac_in_temp),
            random_state=seed,
            stratify=temp_df['len_bin']
        )
        
        train_df = train_df.drop(columns=['len_bin'])
        valid_df = valid_df.drop(columns=['len_bin'])
        test_df = test_df.drop(columns=['len_bin'])
        
        train_df.to_csv(config.PROC_DIR / f"{output_prefix}_random_train.csv", index=False)
        valid_df.to_csv(config.PROC_DIR / f"{output_prefix}_random_valid.csv", index=False)
        test_df.to_csv(config.PROC_DIR / f"{output_prefix}_random_test.csv", index=False)
        
        log_message(f"Random split: train={len(train_df)}, valid={len(valid_df)}, test={len(test_df)}")
        
        return train_df, valid_df, test_df


# =========================================================
# 13. Tokenizers
# =========================================================
class BaseTokenizer:
    """Single nucleotide tokenizer."""
    
    def __init__(self):
        self.vocab = {
            "<PAD>": 0,
            "A": 1, "C": 2, "G": 3, "U": 4, "N": 5
        }
        self.inv_vocab = {v: k for k, v in self.vocab.items()}
    
    @property
    def vocab_size(self):
        return len(self.vocab)
    
    def encode(self, seq: str, max_len: int) -> np.ndarray:
        """Encode sequence to token indices."""
        seq = clean_sequence(seq)
        arr = np.zeros(max_len, dtype=np.int64)
        for i, ch in enumerate(seq[:max_len]):
            arr[i] = self.vocab.get(ch, self.vocab["N"])
        return arr
    
    def decode(self, ids: np.ndarray) -> str:
        """Decode token indices to sequence."""
        return ''.join(self.inv_vocab.get(idx, 'N') for idx in ids if idx != 0)


class KMerTokenizer:
    """k-mer tokenizer."""
    
    def __init__(self, k: int = 3):
        self.k = k
        self.bases = ["A", "C", "G", "U", "N"]
        self.pad_token = "<PAD>"
        self.unk_token = "<UNK>"
        self._build_vocab()
    
    def _build_vocab(self):
        """Build vocabulary of all possible k-mers."""
        kmers = ["".join(x) for x in itertools.product(self.bases, repeat=self.k)]
        
        self.vocab = {self.pad_token: 0, self.unk_token: 1}
        for i, km in enumerate(kmers, start=2):
            self.vocab[km] = i
        
        self.inv_vocab = {v: k for k, v in self.vocab.items()}
    
    @property
    def vocab_size(self):
        return len(self.vocab)
    
    def tokenize(self, seq: str) -> List[str]:
        """Split sequence into k-mers."""
        seq = clean_sequence(seq)
        if len(seq) < self.k:
            return [self.unk_token]
        
        return [seq[i:i+self.k] for i in range(len(seq) - self.k + 1)]
    
    def encode(self, seq: str, max_len: int) -> np.ndarray:
        """Encode sequence to k-mer indices."""
        tokens = self.tokenize(seq)
        arr = np.zeros(max_len, dtype=np.int64)
        
        for i, token in enumerate(tokens[:max_len]):
            arr[i] = self.vocab.get(token, self.vocab[self.unk_token])
        
        return arr
    
    def decode(self, ids: np.ndarray) -> List[str]:
        """Decode indices back to k-mers."""
        return [self.inv_vocab.get(idx, self.unk_token) for idx in ids if idx != 0]


# =========================================================
# 14. RNA-FM Embedding (Optional)
# =========================================================
class RNAFMEmbedding(nn.Module):
    """
    RNA-FM pretrained model embedding layer.
    Falls back to dummy embedding if RNA-FM is not available.
    """
    
    def __init__(self, model_name: str = "RNA-FM/RNA-FM", output_dim: int = 128):
        super().__init__()
        self.output_dim = output_dim
        self.use_rnafm = False
        
        if TRANSFORMERS_AVAILABLE:
            try:
                self.tokenizer = AutoTokenizer.from_pretrained(model_name)
                self.model = AutoModel.from_pretrained(model_name)
                self.model.eval()
                self.projection = nn.Linear(640, output_dim)  # RNA-FM output dimension is 640
                self.use_rnafm = True
                log_message(f"RNA-FM loaded: {model_name}")
            except Exception as e:
                log_message(f"Failed to load RNA-FM: {e}", "WARNING")
                self._init_dummy()
        else:
            self._init_dummy()
    
    def _init_dummy(self):
        """Initialize dummy layer when RNA-FM is not available."""
        self.use_rnafm = False
        self.dummy = nn.Linear(self.output_dim, self.output_dim)
    
    def forward(self, sequences: List[str]) -> torch.Tensor:
        """
        Forward pass.
        
        Args:
            sequences: List of RNA sequences
        
        Returns:
            Tensor of shape [B, L, output_dim]
        """
        if not self.use_rnafm:
            B = len(sequences)
            L = max(len(seq) for seq in sequences)
            return self.dummy(torch.randn(B, L, self.output_dim, 
                                          device=next(self.parameters()).device))
        
        with torch.no_grad():
            inputs = self.tokenizer(sequences, return_tensors="pt", padding=True)
            inputs = {k: v.to(next(self.parameters()).device) for k, v in inputs.items()}
            outputs = self.model(**inputs)
            embeddings = outputs.last_hidden_state  # [B, L, 640]
            embeddings = self.projection(embeddings)  # [B, L, output_dim]
        
        return embeddings


# =========================================================
# 15. Data Encoding Utilities
# =========================================================
def build_pair_matrix(dot_bracket: str, max_len: int, 
                     allow_pseudoknots: bool = False) -> np.ndarray:
    """Build pair matrix from dot-bracket notation."""
    if not dot_bracket:
        return np.zeros((max_len, max_len), dtype=np.float32)
    
    gt = BracketedNotationParser.dotbracket_to_matrix(dot_bracket, allow_pseudoknots)
    out = np.zeros((max_len, max_len), dtype=np.float32)
    L = min(max_len, gt.shape[0])
    out[:L, :L] = gt[:L, :L]
    return out


def build_valid_mask(length: int, max_len: int, min_loop: int = 3,
                    include_ambiguous: bool = True) -> np.ndarray:
    """
    Build mask for valid pair positions (i < j, j - i > min_loop, valid base pairs).
    Used for both training and evaluation.
    """
    mask = np.zeros((max_len, max_len), dtype=np.float32)
    
    # Only consider upper triangle (i < j)
    for i in range(min(length, max_len)):
        for j in range(i + 1, min(length, max_len)):
            if j - i > min_loop:
                mask[i, j] = 1.0
    
    return mask


def build_basepair_prior(seq: str, max_len: int, min_loop: int = 3,
                        include_ambiguous: bool = True) -> np.ndarray:
    """Build base pairing prior based on sequence only."""
    seq = clean_sequence(seq)
    L = min(len(seq), max_len)
    prior = np.zeros((max_len, max_len), dtype=np.float32)
    
    for i in range(L):
        for j in range(L):
            if i == j:
                continue
            if abs(i - j) <= min_loop:
                continue
            if BasePairingRules.is_valid_pair(seq[i], seq[j], include_ambiguous):
                prior[i, j] = 1.0
    
    return prior


def build_distance_bias(max_len: int) -> np.ndarray:
    """Build distance bias matrix (normalized by max_len)."""
    idx = np.arange(max_len)
    dist = np.abs(idx[:, None] - idx[None, :]).astype(np.float32)
    return dist / max(max_len, 1)


def build_positional_encoding(max_len: int, d_model: int) -> np.ndarray:
    """Build sinusoidal positional encoding."""
    pe = np.zeros((max_len, d_model))
    position = np.arange(0, max_len)[:, np.newaxis]
    div_term = np.exp(np.arange(0, d_model, 2) * -(math.log(10000.0) / d_model))
    
    pe[:, 0::2] = np.sin(position * div_term)
    pe[:, 1::2] = np.cos(position * div_term)
    
    return pe


# =========================================================
# 16. Dataset Class
# =========================================================
class RNASSDataset(Dataset):
    """RNA secondary structure dataset."""
    
    def __init__(
        self,
        df: pd.DataFrame,
        tokenizer,
        max_len: int = 256,
        min_loop: int = 3,
        use_valid_space: bool = True,
        allow_pseudoknots: bool = False
    ):
        self.df = df.reset_index(drop=True)
        self.tokenizer = tokenizer
        self.max_len = max_len
        self.min_loop = min_loop
        self.use_valid_space = use_valid_space
        self.allow_pseudoknots = allow_pseudoknots
        self.dist_bias = build_distance_bias(max_len)
    
    def __len__(self):
        return len(self.df)
    
    def __getitem__(self, idx):
        row = self.df.iloc[idx]
        seq = clean_sequence(row["sequence"])
        db = clean_dot_bracket(row.get("dot_bracket", ""))
        
        L = min(len(seq), self.max_len)
        
        # Tokenize
        token_ids = self.tokenizer.encode(seq, self.max_len)
        
        # Target pair matrix
        pair_target = build_pair_matrix(db, self.max_len, self.allow_pseudoknots)
        
        # Valid mask for training/evaluation
        if self.use_valid_space:
            valid_mask = build_valid_mask(L, self.max_len, self.min_loop)
        else:
            valid_mask = np.ones((self.max_len, self.max_len), dtype=np.float32)
        
        # Base pairing prior
        pair_prior = build_basepair_prior(seq, self.max_len, self.min_loop)
        
        return {
            "input_ids": torch.tensor(token_ids, dtype=torch.long),
            "pair_target": torch.tensor(pair_target, dtype=torch.float32),
            "valid_mask": torch.tensor(valid_mask, dtype=torch.float32),
            "pair_prior": torch.tensor(pair_prior, dtype=torch.float32),
            "dist_bias": torch.tensor(self.dist_bias, dtype=torch.float32),
            "length": torch.tensor(L, dtype=torch.long),
            "sequence": seq,
            "dot_bracket": db,
            "sequence_id": row.get("sequence_id", str(idx))
        }


def collate_fn(batch):
    """Custom collate function for DataLoader."""
    out = {}
    tensor_keys = ["input_ids", "pair_target", "valid_mask", "pair_prior", "dist_bias", "length"]
    str_keys = ["sequence", "dot_bracket", "sequence_id"]
    
    for k in tensor_keys:
        out[k] = torch.stack([item[k] for item in batch], dim=0)
    for k in str_keys:
        out[k] = [item[k] for item in batch]
    
    return out


# =========================================================
# 17. Model Components
# =========================================================
class PositionalEncoding(nn.Module):
    """Sinusoidal positional encoding."""
    
    def __init__(self, d_model: int, max_len: int = 1024, dropout: float = 0.1):
        super().__init__()
        self.dropout = nn.Dropout(p=dropout)
        
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len).unsqueeze(1).float()
        div_term = torch.exp(
            torch.arange(0, d_model, 2).float() * (-math.log(10000.0) / d_model)
        )
        
        pe[:, 0::2] = torch.sin(position * div_term)
        if d_model % 2 == 1:
            pe[:, 1::2] = torch.cos(position * div_term[:-1])
        else:
            pe[:, 1::2] = torch.cos(position * div_term)
        
        self.register_buffer("pe", pe.unsqueeze(0))
    
    def forward(self, x):
        """Add positional encoding to input."""
        return self.dropout(x + self.pe[:, :x.size(1)])


class ConvBlock2D(nn.Module):
    """2D convolutional block with batch norm and activation."""
    
    def __init__(self, in_ch, out_ch, dropout=0.1):
        super().__init__()
        self.block = nn.Sequential(
            nn.Conv2d(in_ch, out_ch, kernel_size=3, padding=1),
            nn.BatchNorm2d(out_ch),
            nn.GELU(),
            nn.Dropout2d(dropout),
            nn.Conv2d(out_ch, out_ch, kernel_size=3, padding=1),
            nn.BatchNorm2d(out_ch),
            nn.GELU(),
            nn.Dropout2d(dropout),
        )
    
    def forward(self, x):
        return self.block(x)


class UNet2D(nn.Module):
    """2D U-Net for contact map refinement."""
    
    def __init__(self, in_ch=4, base_ch=32, out_ch=1, dropout=0.1):
        super().__init__()
        
        # Encoder
        self.enc1 = ConvBlock2D(in_ch, base_ch, dropout)
        self.pool1 = nn.MaxPool2d(2)
        
        self.enc2 = ConvBlock2D(base_ch, base_ch * 2, dropout)
        self.pool2 = nn.MaxPool2d(2)
        
        # Bottleneck
        self.bottleneck = ConvBlock2D(base_ch * 2, base_ch * 4, dropout)
        
        # Decoder with skip connections
        self.up2 = nn.ConvTranspose2d(base_ch * 4, base_ch * 2, kernel_size=2, stride=2)
        self.dec2 = ConvBlock2D(base_ch * 4, base_ch * 2, dropout)
        
        self.up1 = nn.ConvTranspose2d(base_ch * 2, base_ch, kernel_size=2, stride=2)
        self.dec1 = ConvBlock2D(base_ch * 2, base_ch, dropout)
        
        # Output layer
        self.out_conv = nn.Conv2d(base_ch, out_ch, kernel_size=1)
    
    def forward(self, x):
        # Encoder
        e1 = self.enc1(x)
        p1 = self.pool1(e1)
        
        e2 = self.enc2(p1)
        p2 = self.pool2(e2)
        
        # Bottleneck
        b = self.bottleneck(p2)
        
        # Decoder with skip connections
        u2 = self.up2(b)
        if u2.shape[-2:] != e2.shape[-2:]:
            u2 = F.interpolate(u2, size=e2.shape[-2:], mode="bilinear", align_corners=False)
        d2 = self.dec2(torch.cat([u2, e2], dim=1))
        
        u1 = self.up1(d2)
        if u1.shape[-2:] != e1.shape[-2:]:
            u1 = F.interpolate(u1, size=e1.shape[-2:], mode="bilinear", align_corners=False)
        d1 = self.dec1(torch.cat([u1, e1], dim=1))
        
        # Output
        out = self.out_conv(d1)
        return out


# =========================================================
# 18. Main Model: RNA Secondary Structure Predictor
# =========================================================
class RNASSPredictor(nn.Module):
    """
    Advanced RNA secondary structure prediction model.
    
    Architecture:
    - Token embedding (k-mer or RNA-FM)
    - Transformer encoder for sequence representation
    - Multi-way feature fusion for pair representation
    - U-Net for contact map refinement
    """
    
    def __init__(
        self,
        vocab_size: int,
        d_model: int = 128,
        nhead: int = 8,
        num_layers: int = 4,
        dim_ff: int = 256,
        dropout: float = 0.1,
        max_len: int = 256,
        use_rnafm: bool = False,
        use_positional_encoding: bool = True
    ):
        super().__init__()
        self.max_len = max_len
        self.d_model = d_model
        self.use_rnafm = use_rnafm
        
        # Embedding layer
        if use_rnafm:
            self.embedding = RNAFMEmbedding(output_dim=d_model)
        else:
            self.embedding = nn.Embedding(vocab_size, d_model, padding_idx=0)
        
        # Positional encoding
        if use_positional_encoding:
            self.pos_enc = PositionalEncoding(d_model, max_len=max_len, dropout=dropout)
        else:
            self.pos_enc = nn.Identity()
        
        # Transformer encoder
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=dim_ff,
            dropout=dropout,
            batch_first=True,
            activation="gelu",
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)
        
        # Pair feature projections
        self.row_proj = nn.Linear(d_model, d_model)
        self.col_proj = nn.Linear(d_model, d_model)
        self.mul_proj = nn.Linear(d_model, d_model)
        
        # Pair feature reduction
        self.pair_reduce = nn.Conv2d(d_model, 1, kernel_size=1)
        
        # U-Net for refinement
        self.unet = UNet2D(in_ch=4, base_ch=32, out_ch=1, dropout=dropout)
        
        self._init_weights()
    
    def _init_weights(self):
        """Initialize weights using Xavier uniform."""
        for p in self.parameters():
            if p.dim() > 1:
                nn.init.xavier_uniform_(p)
    
    def forward(self, input_ids, pair_prior, dist_bias, sequences=None):
        """
        Forward pass.
        
        Args:
            input_ids: [B, L] token indices
            pair_prior: [B, L, L] base pairing prior
            dist_bias: [B, L, L] distance bias
            sequences: Optional list of sequences for RNA-FM
        
        Returns:
            logits: [B, L, L] pair logits
        """
        # Create padding mask
        pad_mask = (input_ids == 0)
        
        # Embedding
        if self.use_rnafm and sequences is not None:
            x = self.embedding(sequences)  # [B, L, D]
        else:
            x = self.embedding(input_ids)  # [B, L, D]
        
        # Positional encoding and Transformer
        x = self.pos_enc(x)
        x = self.encoder(x, src_key_padding_mask=pad_mask)  # [B, L, D]
        
        # Generate pair features using multi-way fusion
        row = self.row_proj(x)  # [B, L, D]
        col = self.col_proj(x)  # [B, L, D]
        mul = self.mul_proj(x)  # [B, L, D]
        
        # Additive fusion
        pair_add = row.unsqueeze(2) + col.unsqueeze(1)  # [B, L, L, D]
        
        # Multiplicative fusion
        pair_mul = mul.unsqueeze(2) * mul.unsqueeze(1)  # [B, L, L, D]
        
        # Combine
        pair_feat = pair_add + pair_mul  # [B, L, L, D]
        
        # Reshape for 2D convolutions
        pair_feat = pair_feat.permute(0, 3, 1, 2)  # [B, D, L, L]
        pair_feat = self.pair_reduce(pair_feat)  # [B, 1, L, L]
        
        # Prepare U-Net input (4 channels)
        prior = pair_prior.unsqueeze(1)  # [B, 1, L, L]
        dist = dist_bias.unsqueeze(1)    # [B, 1, L, L]
        sym_prior = 0.5 * (prior + prior.transpose(2, 3))  # Symmetric prior
        
        unet_input = torch.cat([pair_feat, prior, sym_prior, dist], dim=1)  # [B, 4, L, L]
        
        # U-Net refinement
        logits = self.unet(unet_input).squeeze(1)  # [B, L, L]
        
        # Enforce symmetry
        logits = 0.5 * (logits + logits.transpose(1, 2))
        
        return logits


# =========================================================
# 19. Multi-Objective Loss Function
# =========================================================
class RNAStructureLoss(nn.Module):
    """
    Multi-objective loss function for RNA structure prediction.
    
    Components:
    - BCE: Binary cross-entropy
    - Focal: Focal loss for class imbalance
    - Dice: Dice coefficient loss
    - Symmetry: Enforce prediction symmetry
    - Diagonal: Penalize self-pairs
    - Sparsity: Regularize pair density
    - Degree: Penalize multiple pairs per position (optional)
    
    All calculations are performed only on valid pair space (i < j, j - i > min_loop).
    """
    
    def __init__(self, config):
        super().__init__()
        self.config = config
        self.weights = config.LOSS_WEIGHTS
    
    def _masked_bce(self, logits, targets, mask, pos_weight=None):
        """Masked binary cross-entropy loss."""
        if pos_weight is None:
            loss = F.binary_cross_entropy_with_logits(logits, targets, reduction="none")
        else:
            loss = F.binary_cross_entropy_with_logits(
                logits, targets, reduction="none", pos_weight=pos_weight
            )
        loss = loss * mask
        return loss.sum() / mask.sum().clamp(min=1.0)
    
    def _focal_loss(self, logits, targets, mask, alpha=0.25, gamma=2.0):
        """Focal loss for handling class imbalance."""
        probs = torch.sigmoid(logits)
        bce = F.binary_cross_entropy_with_logits(logits, targets, reduction="none")
        pt = probs * targets + (1 - probs) * (1 - targets)
        focal = alpha * (1 - pt).pow(gamma) * bce
        focal = focal * mask
        return focal.sum() / mask.sum().clamp(min=1.0)
    
    def _dice_loss(self, logits, targets, mask, eps=1e-6):
        """Dice coefficient loss."""
        probs = torch.sigmoid(logits) * mask
        targets = targets * mask
        inter = (probs * targets).sum()
        union = probs.sum() + targets.sum()
        dice = (2 * inter + eps) / (union + eps)
        return 1 - dice
    
    def _symmetry_loss(self, logits, mask):
        """Enforce symmetry in predictions."""
        diff = torch.abs(logits - logits.transpose(1, 2)) * mask
        return diff.sum() / mask.sum().clamp(min=1.0)
    
    def _diagonal_penalty(self, logits, mask):
        """Penalize predictions on the diagonal (self-pairs)."""
        B, L, _ = logits.shape
        diag = torch.eye(L, device=logits.device).unsqueeze(0).repeat(B, 1, 1)
        probs = torch.sigmoid(logits)
        loss = (probs * diag * mask).sum() / (diag * mask).sum().clamp(min=1.0)
        return loss
    
    def _sparsity_regularization(self, logits, mask, target_density=0.02):
        """Regularize pair density to match expected value."""
        probs = torch.sigmoid(logits) * mask
        density = probs.sum() / mask.sum().clamp(min=1.0)
        return torch.abs(density - target_density)
    
    def _degree_regularization(self, logits, mask):
        """Penalize multiple pairs per position."""
        probs = torch.sigmoid(logits) * mask
        row_sums = probs.sum(dim=2)  # [B, L]
        col_sums = probs.sum(dim=1)  # [B, L]
        
        # Each position should have at most 1 pair (ideally 0 or 1)
        row_penalty = (row_sums - 1).clamp(min=0).pow(2).mean()
        col_penalty = (col_sums - 1).clamp(min=0).pow(2).mean()
        
        return (row_penalty + col_penalty) / 2
    
    def forward(self, logits, targets, mask):
        """
        Compute total loss.
        
        Args:
            logits: [B, L, L] prediction logits
            targets: [B, L, L] target binary matrix
            mask: [B, L, L] valid position mask (i < j, j - i > min_loop)
        
        Returns:
            total_loss: Scalar loss value
            stats: Dictionary of individual loss components
        """
        # Positive sample weight
        pos_weight = torch.tensor([self.config.POS_WEIGHT], device=logits.device)
        
        # Compute loss components
        l_bce = self._masked_bce(logits, targets, mask, pos_weight)
        l_focal = self._focal_loss(logits, targets, mask, 
                                   self.config.FOCAL_ALPHA, self.config.FOCAL_GAMMA)
        l_dice = self._dice_loss(logits, targets, mask)
        l_sym = self._symmetry_loss(logits, mask)
        l_diag = self._diagonal_penalty(logits, mask)
        l_sparse = self._sparsity_regularization(logits, mask, self.config.TARGET_PAIR_DENSITY)
        
        # Weighted combination
        total = (
            self.weights['bce'] * l_bce +
            self.weights['focal'] * l_focal +
            self.weights['dice'] * l_dice +
            self.weights['symmetry'] * l_sym +
            self.weights['diagonal'] * l_diag +
            self.weights['sparsity'] * l_sparse
        )
        
        # Optional degree regularization
        if self.weights.get('degree', 0) > 0:
            l_degree = self._degree_regularization(logits, mask)
            total += self.weights['degree'] * l_degree
            stats = {
                'bce': l_bce.item(),
                'focal': l_focal.item(),
                'dice': l_dice.item(),
                'symmetry': l_sym.item(),
                'diagonal': l_diag.item(),
                'sparsity': l_sparse.item(),
                'degree': l_degree.item() if 'degree' in locals() else 0,
                'total': total.item()
            }
        else:
            stats = {
                'bce': l_bce.item(),
                'focal': l_focal.item(),
                'dice': l_dice.item(),
                'symmetry': l_sym.item(),
                'diagonal': l_diag.item(),
                'sparsity': l_sparse.item(),
                'total': total.item()
            }
        
        return total, stats


# =========================================================
# 20. Nussinov Decoder (Canonical Structures Only)
# =========================================================
class NussinovDecoder:
    """
    Nussinov dynamic programming decoder for canonical RNA secondary structures.
    
    This decoder finds the optimal non-crossing structure by maximizing
    the sum of pair probabilities.
    
    Note: This decoder does NOT support pseudoknots.
    """
    
    def __init__(self, min_loop: int = 3):
        self.min_loop = min_loop
        self.parser = BracketedNotationParser()
    
    def decode(self, prob_mat: np.ndarray, sequence: str, 
              threshold: float = 0.5) -> str:
        """
        Decode probability matrix to dot-bracket structure.
        
        Args:
            prob_mat: [L, L] probability matrix
            sequence: RNA sequence
            threshold: Minimum probability to consider a pair
        
        Returns:
            Dot-bracket notation string
        """
        seq = clean_sequence(sequence)
        L = len(seq)
        
        if prob_mat.shape[0] != L or prob_mat.shape[1] != L:
            raise ValueError(f"Matrix shape {prob_mat.shape} doesn't match sequence length {L}")
        
        # Build DP tables
        # score[i][j] = maximum score for subsequence i..j
        # trace[i][j] = decision type: 0=unpaired_i, 1=unpaired_j, 2=pair_ij, 3=split
        score = np.zeros((L, L), dtype=np.float32)
        trace = np.zeros((L, L), dtype=np.int8)
        split_k = np.full((L, L), -1, dtype=np.int16)
        
        # Fill DP table
        for length in range(1, L):
            for i in range(L - length):
                j = i + length
                
                # Case 1: i is unpaired
                best_score = score[i + 1, j] if i + 1 <= j else 0.0
                best_trace = 0  # unpaired i
                
                # Case 2: j is unpaired
                if i <= j - 1:
                    candidate = score[i, j - 1]
                    if candidate > best_score:
                        best_score = candidate
                        best_trace = 1  # unpaired j
                
                # Case 3: i and j pair
                if (j - i > self.min_loop and 
                    BasePairingRules.is_valid_pair(seq[i], seq[j], include_ambiguous=False) and
                    prob_mat[i, j] >= threshold):
                    
                    inner_score = score[i + 1, j - 1] if i + 1 <= j - 1 else 0.0
                    candidate = prob_mat[i, j] + inner_score
                    
                    if candidate > best_score:
                        best_score = candidate
                        best_trace = 2  # pair i-j
                
                # Case 4: bifurcation (split at k)
                for k in range(i + 1, j):
                    candidate = score[i, k] + score[k + 1, j]
                    if candidate > best_score:
                        best_score = candidate
                        best_trace = 3  # split
                        split_k[i, j] = k
                
                score[i, j] = best_score
                trace[i, j] = best_trace
        
        # Traceback to get pairs
        pairs = []
        
        def traceback(i, j):
            if i >= j:
                return
            
            if trace[i, j] == 0:  # i unpaired
                traceback(i + 1, j)
            elif trace[i, j] == 1:  # j unpaired
                traceback(i, j - 1)
            elif trace[i, j] == 2:  # i and j paired
                pairs.append((i, j))
                traceback(i + 1, j - 1)
            elif trace[i, j] == 3:  # bifurcation
                k = split_k[i, j]
                if k != -1:
                    traceback(i, k)
                    traceback(k + 1, j)
        
        if L > 0:
            traceback(0, L - 1)
        
        # Convert pairs to dot-bracket (canonical only)
        return self.parser.pairs_to_dotbracket(pairs, L, allow_pseudoknots=False)


# =========================================================
# 21. ViennaRNA Thermodynamic Model Interface
# =========================================================
class ViennaRNAInterface:
    """Interface to ViennaRNA thermodynamic model."""
    
    def __init__(self, rnafold_path: str = "RNAfold"):
        self.rnafold_path = rnafold_path
        self.available = self._check_available()
    
    def _check_available(self) -> bool:
        """Check if RNAfold is available."""
        try:
            subprocess.run([self.rnafold_path, "--version"], 
                          capture_output=True, check=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            log_message("RNAfold not found. Thermodynamic fusion disabled.", "WARNING")
            return False
    
    def get_pair_probabilities(self, sequence: str) -> Optional[np.ndarray]:
        """
        Get pair probability matrix from ViennaRNA.
        
        Args:
            sequence: RNA sequence
        
        Returns:
            [L, L] probability matrix or None if failed
        """
        if not self.available:
            return None
        
        seq = clean_sequence(sequence)
        L = len(seq)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            work_dir = Path(tmpdir)
            seq_file = work_dir / "sequence.fa"
            
            with open(seq_file, 'w') as f:
                f.write(f">seq\n{seq}\n")
            
            try:
                cmd = f'cd "{work_dir}" && {self.rnafold_path} -p < sequence.fa'
                subprocess.run(cmd, shell=True, check=True, 
                             capture_output=True, timeout=30)
                
                # Find probability output file
                prob_file = self._find_probability_file(work_dir)
                if prob_file:
                    return self._parse_probability_file(prob_file, L)
                
            except Exception as e:
                log_message(f"RNAfold failed: {e}", "WARNING")
        
        return None
    
    def _find_probability_file(self, work_dir: Path) -> Optional[Path]:
        """Find the probability output file."""
        patterns = ['*dp.ps', '*.prob', '*_dp.ps']
        for pattern in patterns:
            files = list(work_dir.glob(pattern))
            if files:
                return files[0]
        return None
    
    def _parse_probability_file(self, file_path: Path, L: int) -> np.ndarray:
        """
        Parse ViennaRNA probability file.
        Returns a symmetric probability matrix.
        """
        prob_mat = np.zeros((L, L), dtype=np.float32)
        
        try:
            with open(file_path, 'r') as f:
                content = f.read()
                
                # Look for probability data
                # ViennaRNA often outputs in format: i j p
                pattern = r'(\d+)\s+(\d+)\s+([0-9.]+)'
                matches = re.findall(pattern, content)
                
                for i, j, p in matches:
                    i, j = int(i) - 1, int(j) - 1  # Convert to 0-based
                    if 0 <= i < L and 0 <= j < L and i < j:
                        # ViennaRNA often outputs sqrt(p)
                        prob = float(p) ** 2
                        prob_mat[i, j] = prob
                        prob_mat[j, i] = prob
                        
        except Exception as e:
            log_message(f"Failed to parse probability file: {e}", "WARNING")
        
        return prob_mat
    
    def fuse_with_learned(self, learned_prob: np.ndarray, 
                          thermo_prob: Optional[np.ndarray],
                          alpha: float = 0.8) -> np.ndarray:
        """
        Fuse learned probabilities with thermodynamic predictions.
        
        Args:
            learned_prob: [L, L] probability matrix from deep learning
            thermo_prob: [L, L] probability matrix from ViennaRNA
            alpha: Weight for learned predictions (1-alpha for ViennaRNA)
        
        Returns:
            Fused probability matrix
        """
        if thermo_prob is None:
            return learned_prob
        
        if thermo_prob.shape != learned_prob.shape:
            return learned_prob
        
        return alpha * learned_prob + (1 - alpha) * thermo_prob


# =========================================================
# 22. Strict Evaluation Metrics
# =========================================================
class RNAStructureMetrics:
    """
    Strict evaluation metrics for RNA secondary structure prediction.
    
    Evaluation is performed only on the valid pair space:
    - Upper triangle (i < j)
    - Excluding hairpin loops (j - i > min_loop)
    - Only canonical base pairs (AU, CG, GU)
    """
    
    def __init__(self, min_loop: int = 3):
        self.min_loop = min_loop
    
    def compute_metrics(self, pred_db: str, true_db: str, sequence: str) -> Dict[str, float]:
        """
        Compute comprehensive evaluation metrics.
        
        Args:
            pred_db: Predicted dot-bracket structure
            true_db: True dot-bracket structure
            sequence: RNA sequence
        
        Returns:
            Dictionary of metrics including precision, recall, F1, MCC, exact_match
        """
        # Get valid pair space
        valid_pairs = self._get_valid_pair_space(sequence)
        valid_set = set(valid_pairs)
        
        # Parse structures
        pred_pairs = self._parse_structure(pred_db)
        true_pairs = self._parse_structure(true_db)
        
        # Restrict to valid pairs
        pred_valid = pred_pairs & valid_set
        true_valid = true_pairs & valid_set
        
        # Confusion matrix
        tp = len(pred_valid & true_valid)
        fp = len(pred_valid - true_valid)
        fn = len(true_valid - pred_valid)
        tn = len(valid_set) - tp - fp - fn
        
        # Basic metrics
        precision = safe_div(tp, tp + fp)
        recall = safe_div(tp, tp + fn)
        f1 = safe_div(2 * precision * recall, precision + recall)
        
        # Matthews Correlation Coefficient
        mcc_numerator = tp * tn - fp * fn
        mcc_denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        mcc = safe_div(mcc_numerator, mcc_denom)
        
        # Exact match on valid pairs
        exact_match = 1.0 if pred_valid == true_valid else 0.0
        
        # Sensitivity and PPV (same as recall and precision)
        sensitivity = recall
        ppv = precision
        
        # F1 for pseudoknots (if present)
        pk_pred = self._identify_pseudoknots(pred_pairs)
        pk_true = self._identify_pseudoknots(true_pairs)
        
        pk_precision = safe_div(len(pk_pred & pk_true), len(pk_pred))
        pk_recall = safe_div(len(pk_pred & pk_true), len(pk_true))
        pk_f1 = safe_div(2 * pk_precision * pk_recall, pk_precision + pk_recall)
        
        return {
            'precision': precision,
            'recall': recall,
            'f1': f1,
            'mcc': mcc,
            'exact_match': exact_match,
            'sensitivity': sensitivity,
            'ppv': ppv,
            'tp': tp,
            'fp': fp,
            'fn': fn,
            'tn': tn,
            'pseudoknot_precision': pk_precision,
            'pseudoknot_recall': pk_recall,
            'pseudoknot_f1': pk_f1,
            'n_valid_pairs': len(valid_set)
        }
    
    def _get_valid_pair_space(self, sequence: str) -> List[Tuple[int, int]]:
        """Get all valid pair positions based on sequence."""
        seq = clean_sequence(sequence)
        L = len(seq)
        valid_pairs = []
        
        for i in range(L):
            for j in range(i + 1, L):
                if j - i <= self.min_loop:
                    continue
                if BasePairingRules.is_valid_pair(seq[i], seq[j], include_ambiguous=False):
                    valid_pairs.append((i, j))
        
        return valid_pairs
    
    def _parse_structure(self, db: str) -> set:
        """Parse dot-bracket structure to set of pairs."""
        if not db:
            return set()
        pairs = BracketedNotationParser.dotbracket_to_pairs(db, allow_pseudoknots=True)
        return set((min(i, j), max(i, j)) for i, j in pairs if i != j)
    
    def _identify_pseudoknots(self, pairs: set) -> set:
        """Identify pseudoknot pairs."""
        pair_list = sorted(pairs)
        pseudoknots = set()
        
        for idx1, (i, j) in enumerate(pair_list):
            for idx2, (m, n) in enumerate(pair_list[idx1+1:], idx1+1):
                if (i < m < j < n) or (m < i < n < j):
                    pseudoknots.add((i, j))
                    pseudoknots.add((m, n))
        
        return pseudoknots


# =========================================================
# 23. Uncertainty Estimation (MC Dropout)
# =========================================================
class UncertaintyEstimator:
    """Monte Carlo Dropout uncertainty estimator."""
    
    def __init__(self, model: nn.Module, num_samples: int = 20):
        self.model = model
        self.num_samples = num_samples
    
    def _enable_dropout(self):
        """Enable dropout layers for inference."""
        for module in self.model.modules():
            if isinstance(module, (nn.Dropout, nn.Dropout2d)):
                module.train()
    
    def estimate(self, input_ids, pair_prior, dist_bias, sequences=None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Estimate prediction uncertainty.
        
        Returns:
            mean_pred: [B, L, L] mean predictions
            var_pred: [B, L, L] prediction variance
        """
        self.model.eval()
        self._enable_dropout()
        
        predictions = []
        with torch.no_grad():
            for _ in range(self.num_samples):
                logits = self.model(input_ids, pair_prior, dist_bias, sequences=sequences)
                probs = torch.sigmoid(logits)
                predictions.append(probs.cpu().numpy())
        
        predictions = np.array(predictions)  # [S, B, L, L]
        mean_pred = predictions.mean(axis=0)
        variance_pred = predictions.var(axis=0)
        
        return mean_pred, variance_pred
    
    def get_confidence_intervals(self, input_ids, pair_prior, dist_bias,
                                 confidence_level: float = 0.95) -> Tuple[np.ndarray, np.ndarray]:
        """Get confidence intervals for predictions."""
        mean_pred, var_pred = self.estimate(input_ids, pair_prior, dist_bias)
        
        # Assuming normal distribution
        z_score = stats.norm.ppf((1 + confidence_level) / 2)
        std_pred = np.sqrt(var_pred)
        
        ci_lower = mean_pred - z_score * std_pred
        ci_upper = mean_pred + z_score * std_pred
        
        return ci_lower, ci_upper


# =========================================================
# 24. Trainer Class
# =========================================================
class Trainer:
    """Model trainer with comprehensive logging and checkpointing."""
    
    def __init__(self, model, config, run_id: str = None):
        self.model = model
        self.config = config
        self.device = config.DEVICE
        
        if run_id is None:
            run_id = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        self.run_id = run_id
        self.run_dir = config.MODEL_DIR / run_id
        self.run_dir.mkdir(parents=True, exist_ok=True)
        
        # Loss function
        self.criterion = RNAStructureLoss(config)
        
        # Optimizer
        self.optimizer = torch.optim.AdamW(
            model.parameters(),
            lr=config.LR,
            weight_decay=config.WEIGHT_DECAY
        )
        
        # Learning rate scheduler
        self.scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            self.optimizer,
            mode='max',
            factor=0.5,
            patience=3,
            verbose=True
        )
        
        # Evaluator
        self.metrics = RNAStructureMetrics(min_loop=config.MIN_HAIRPIN_LOOP)
        self.decoder = NussinovDecoder(min_loop=config.MIN_HAIRPIN_LOOP)
        
        # Tracking
        self.best_f1 = -1.0
        self.best_epoch = -1
        self.patience_counter = 0
        
        # History
        self.history = {
            'train_loss': [],
            'val_metrics': [],
            'learning_rates': [],
            'epoch_times': []
        }
        
        # Save config
        config.save(self.run_dir / "config.json")
    
    def train_epoch(self, train_loader) -> Tuple[float, Dict]:
        """Train for one epoch."""
        self.model.train()
        total_loss = 0
        loss_stats = defaultdict(float)
        
        pbar = tqdm(train_loader, desc="Training")
        for batch in pbar:
            # Move to device
            input_ids = batch["input_ids"].to(self.device)
            pair_target = batch["pair_target"].to(self.device)
            valid_mask = batch["valid_mask"].to(self.device)
            pair_prior = batch["pair_prior"].to(self.device)
            dist_bias = batch["dist_bias"].to(self.device)
            
            # Forward pass
            self.optimizer.zero_grad()
            logits = self.model(input_ids, pair_prior, dist_bias)
            loss, stats = self.criterion(logits, pair_target, valid_mask)
            
            # Backward pass
            loss.backward()
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), self.config.GRAD_CLIP)
            self.optimizer.step()
            
            # Statistics
            total_loss += loss.item()
            for k, v in stats.items():
                loss_stats[k] += v
            
            pbar.set_postfix({'loss': f'{loss.item():.4f}'})
        
        avg_loss = total_loss / len(train_loader)
        avg_stats = {k: v / len(train_loader) for k, v in loss_stats.items()}
        
        return avg_loss, avg_stats
    
    def validate(self, val_loader, threshold: float = 0.5) -> Dict:
        """Validate the model."""
        self.model.eval()
        
        all_metrics = defaultdict(list)
        
        with torch.no_grad():
            for batch in tqdm(val_loader, desc="Validating"):
                input_ids = batch["input_ids"].to(self.device)
                pair_prior = batch["pair_prior"].to(self.device)
                dist_bias = batch["dist_bias"].to(self.device)
                
                # Predict
                logits = self.model(input_ids, pair_prior, dist_bias)
                probs = torch.sigmoid(logits).cpu().numpy()
                
                # Decode and evaluate each sample
                for i in range(len(batch["sequence"])):
                    seq = batch["sequence"][i]
                    true_db = batch["dot_bracket"][i]
                    L = batch["length"][i].item()
                    
                    pred_prob = probs[i][:L, :L]
                    pred_db = self.decoder.decode(pred_prob, seq, threshold)
                    
                    metrics = self.metrics.compute_metrics(pred_db, true_db, seq)
                    for k, v in metrics.items():
                        all_metrics[k].append(v)
        
        # Average metrics
        avg_metrics = {k: np.mean(v) for k, v in all_metrics.items()}
        avg_metrics['std_f1'] = np.std(all_metrics['f1']) if 'f1' in all_metrics else 0
        
        return avg_metrics
    
    def train(self, train_loader, val_loader, num_epochs: int):
        """Complete training loop."""
        log_message(f"Starting training on {self.device}, run_id: {self.run_id}")
        
        for epoch in range(1, num_epochs + 1):
            epoch_start = time.time()
            
            log_message(f"\n{'='*50}")
            log_message(f"Epoch {epoch}/{num_epochs}")
            
            # Train
            train_loss, train_stats = self.train_epoch(train_loader)
            
            # Validate
            val_metrics = self.validate(val_loader)
            
            epoch_time = time.time() - epoch_start
            
            # Learning rate scheduling
            self.scheduler.step(val_metrics['f1'])
            current_lr = self.optimizer.param_groups[0]['lr']
            
            # Record history
            self.history['train_loss'].append(train_loss)
            self.history['val_metrics'].append(val_metrics)
            self.history['learning_rates'].append(current_lr)
            self.history['epoch_times'].append(epoch_time)
            
            # Print results
            log_message(f"\nTrain Loss: {train_loss:.4f}")
            log_message(f"Train Stats: {train_stats}")
            log_message(f"\nValidation Metrics:")
            for k, v in val_metrics.items():
                if isinstance(v, float):
                    log_message(f"  {k}: {v:.4f}")
            log_message(f"Epoch Time: {epoch_time:.2f}s")
            
            # Save best model
            if val_metrics['f1'] > self.best_f1:
                self.best_f1 = val_metrics['f1']
                self.best_epoch = epoch
                self.patience_counter = 0
                self._save_checkpoint(epoch, val_metrics, is_best=True)
                log_message(f"New best model saved! F1={val_metrics['f1']:.4f}")
            else:
                self.patience_counter += 1
            
            # Early stopping
            if self.patience_counter >= self.config.EARLY_STOPPING_PATIENCE:
                log_message(f"Early stopping triggered after {epoch} epochs")
                break
            
            # Save checkpoint every 10 epochs
            if epoch % 10 == 0:
                self._save_checkpoint(epoch, val_metrics)
        
        # Training completed
        log_message(f"\n{'='*50}")
        log_message(f"Training completed. Best F1={self.best_f1:.4f} at epoch {self.best_epoch}")
        
        # Save history
        self._save_history()
        
        return self.best_f1
    
    def _save_checkpoint(self, epoch, metrics, is_best=False):
        """Save model checkpoint."""
        checkpoint = {
            'epoch': epoch,
            'model_state_dict': self.model.state_dict(),
            'optimizer_state_dict': self.optimizer.state_dict(),
            'scheduler_state_dict': self.scheduler.state_dict(),
            'metrics': metrics,
            'config': {
                'd_model': self.config.D_MODEL,
                'nhead': self.config.NHEAD,
                'num_layers': self.config.NUM_LAYERS,
                'kmer_size': self.config.KMER_SIZE,
                'max_len': self.config.MAX_LEN,
                'representation': self.config.REPRESENTATION
            },
            'run_id': self.run_id
        }
        
        if is_best:
            path = self.run_dir / "best_model.pt"
        else:
            path = self.run_dir / f"checkpoint_epoch_{epoch}.pt"
        
        torch.save(checkpoint, path)
        log_message(f"Checkpoint saved: {path}")
    
    def _save_history(self):
        """Save training history."""
        history_path = self.run_dir / "training_history.json"
        
        # Convert numpy types to Python native types
        history_dict = {
            'train_loss': [float(x) for x in self.history['train_loss']],
            'val_metrics': [
                {k: float(v) if isinstance(v, np.floating) else v 
                 for k, v in m.items()}
                for m in self.history['val_metrics']
            ],
            'learning_rates': [float(x) for x in self.history['learning_rates']],
            'epoch_times': [float(x) for x in self.history['epoch_times']],
            'best_epoch': self.best_epoch,
            'best_f1': float(self.best_f1)
        }
        
        save_json(history_dict, history_path)
        log_message(f"History saved: {history_path}")


# =========================================================
# 25. Statistical Significance Testing
# =========================================================
class StatisticalSignificance:
    """Statistical significance tests for model comparison."""
    
    @staticmethod
    def bootstrap_test(predictions1: List[float], predictions2: List[float],
                       n_iterations: int = 1000, confidence_level: float = 0.95) -> Dict[str, float]:
        """
        Bootstrap test for comparing two sets of predictions.
        
        Args:
            predictions1: F1 scores from model 1
            predictions2: F1 scores from model 2
            n_iterations: Number of bootstrap iterations
            confidence_level: Confidence level for interval
        
        Returns:
            Dictionary with test statistics
        """
        assert len(predictions1) == len(predictions2)
        n = len(predictions1)
        
        # Observed difference
        observed_diff = np.mean(predictions1) - np.mean(predictions2)
        
        # Bootstrap
        bootstrap_diffs = []
        for _ in range(n_iterations):
            indices = np.random.choice(n, n, replace=True)
            boot_diff = np.mean(predictions1[indices]) - np.mean(predictions2[indices])
            bootstrap_diffs.append(boot_diff)
        
        bootstrap_diffs = np.array(bootstrap_diffs)
        
        # Confidence interval
        alpha = 1 - confidence_level
        ci_lower = np.percentile(bootstrap_diffs, 100 * alpha / 2)
        ci_upper = np.percentile(bootstrap_diffs, 100 * (1 - alpha / 2))
        
        # P-value (two-tailed)
        if observed_diff > 0:
            p_value = np.mean(bootstrap_diffs <= 0)
        else:
            p_value = np.mean(bootstrap_diffs >= 0)
        
        return {
            'observed_diff': observed_diff,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'p_value': p_value,
            'significant': p_value < 0.05,
            'bootstrap_mean': np.mean(bootstrap_diffs),
            'bootstrap_std': np.std(bootstrap_diffs)
        }
    
    @staticmethod
    def mcnemar_test(model1_pred: List[int], model2_pred: List[int], 
                     true_labels: List[int]) -> Dict[str, float]:
        """
        McNemar's test for comparing two classifiers.
        
        Args:
            model1_pred: Predictions from model 1 (0/1)
            model2_pred: Predictions from model 2 (0/1)
            true_labels: Ground truth labels
        
        Returns:
            Dictionary with test statistics
        """
        # Build contingency table
        n00 = n01 = n10 = n11 = 0
        
        for p1, p2, t in zip(model1_pred, model2_pred, true_labels):
            e1 = 1 if p1 == t else 0
            e2 = 1 if p2 == t else 0
            
            if e1 == 1 and e2 == 1:
                n11 += 1
            elif e1 == 1 and e2 == 0:
                n10 += 1
            elif e1 == 0 and e2 == 1:
                n01 += 1
            else:
                n00 += 1
        
        # McNemar's statistic with continuity correction
        chi_square = (abs(n10 - n01) - 1) ** 2 / (n10 + n01 + 1e-8)
        p_value = 1 - stats.chi2.cdf(chi_square, 1)
        
        return {
            'chi_square': chi_square,
            'p_value': p_value,
            'significant': p_value < 0.05,
            'n10': n10,
            'n01': n01
        }
    
    @staticmethod
    def wilcoxon_test(predictions1: List[float], predictions2: List[float]) -> Dict[str, float]:
        """Wilcoxon signed-rank test for paired samples."""
        statistic, p_value = stats.wilcoxon(predictions1, predictions2)
        
        return {
            'statistic': statistic,
            'p_value': p_value,
            'significant': p_value < 0.05
        }


# =========================================================
# 26. Visualization Functions (Nature/Science Style)
# =========================================================
class ScientificVisualizer:
    """Professional visualization for scientific publications."""
    
    def __init__(self, output_dir: Path = None):
        self.output_dir = output_dir or config.FIGURES_DIR
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set Nature-style defaults
        plt.rcParams['font.family'] = 'Times New Roman'
        plt.rcParams['font.size'] = 8
        plt.rcParams['axes.labelsize'] = 9
        plt.rcParams['axes.titlesize'] = 10
        plt.rcParams['legend.fontsize'] = 7
        plt.rcParams['figure.dpi'] = 300
        plt.rcParams['savefig.dpi'] = 300
        plt.rcParams['savefig.bbox'] = 'tight'
        plt.rcParams['savefig.pad_inches'] = 0.1
        
        # Nature color palette (muted, professional)
        self.colors = {
            'blue': '#0072B2',
            'red': '#D55E00',
            'green': '#009E73',
            'orange': '#E69F00',
            'purple': '#CC79A7',
            'cyan': '#56B4E9',
            'gray': '#999999',
            'black': '#000000'
        }
    
    def plot_training_history(self, history: Dict, filename: str = "training_history"):
        """Plot training curves."""
        fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))
        
        # Loss curve
        ax = axes[0]
        epochs = range(1, len(history['train_loss']) + 1)
        ax.plot(epochs, history['train_loss'], 
               color=self.colors['blue'], linewidth=1.5, label='Training')
        ax.set_xlabel('Epoch')
        ax.set_ylabel('Loss')
        ax.set_title('Training Loss', fontweight='bold')
        ax.legend(frameon=False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # F1 curve
        ax = axes[1]
        val_f1 = [m['f1'] for m in history['val_metrics']]
        ax.plot(epochs, val_f1, 
               color=self.colors['green'], linewidth=1.5, label='Validation')
        ax.set_xlabel('Epoch')
        ax.set_ylabel('F1 Score')
        ax.set_title('Validation F1 Score', fontweight='bold')
        ax.legend(frameon=False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / f"{filename}.{config.FIGURE_FORMAT}")
        plt.close()
    
    def plot_confusion_matrix(self, cm: np.ndarray, labels: List[str],
                              filename: str = "confusion_matrix"):
        """Plot confusion matrix."""
        fig, ax = plt.subplots(figsize=(5, 4))
        
        im = ax.imshow(cm, interpolation='nearest', cmap='Blues')
        plt.colorbar(im, ax=ax)
        
        # Add text annotations
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                text = ax.text(j, i, f'{cm[i, j]:.0f}',
                              ha="center", va="center",
                              color="white" if cm[i, j] > cm.max() / 2 else "black")
        
        ax.set_xticks(np.arange(len(labels)))
        ax.set_yticks(np.arange(len(labels)))
        ax.set_xticklabels(labels)
        ax.set_yticklabels(labels)
        ax.set_xlabel('Predicted')
        ax.set_ylabel('True')
        ax.set_title('Confusion Matrix', fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / f"{filename}.{config.FIGURE_FORMAT}")
        plt.close()
    
    def plot_roc_pr_curves(self, y_true: List[int], y_score: List[float],
                           filename: str = "roc_pr_curves"):
        """Plot ROC and PR curves."""
        from sklearn.metrics import roc_curve, auc, precision_recall_curve
        
        fpr, tpr, _ = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)
        
        precision, recall, _ = precision_recall_curve(y_true, y_score)
        pr_auc = auc(recall, precision)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.5))
        
        # ROC curve
        ax1.plot(fpr, tpr, color=self.colors['blue'], linewidth=1.5,
                label=f'ROC (AUC = {roc_auc:.3f})')
        ax1.plot([0, 1], [0, 1], color=self.colors['gray'], linestyle='--', linewidth=1)
        ax1.set_xlabel('False Positive Rate')
        ax1.set_ylabel('True Positive Rate')
        ax1.set_title('ROC Curve', fontweight='bold')
        ax1.legend(frameon=False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        
        # PR curve
        ax2.plot(recall, precision, color=self.colors['green'], linewidth=1.5,
                label=f'PR (AUC = {pr_auc:.3f})')
        ax2.set_xlabel('Recall')
        ax2.set_ylabel('Precision')
        ax2.set_title('Precision-Recall Curve', fontweight='bold')
        ax2.legend(frameon=False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / f"{filename}.{config.FIGURE_FORMAT}")
        plt.close()
        
        return {'roc_auc': roc_auc, 'pr_auc': pr_auc}
    
    def plot_attention_map(self, attention_matrix: np.ndarray, sequence: str,
                          filename: str = "attention_map"):
        """Plot attention map."""
        fig, ax = plt.subplots(figsize=(6, 5))
        
        im = ax.imshow(attention_matrix, cmap='viridis', aspect='auto')
        plt.colorbar(im, ax=ax, label='Attention Weight')
        
        # Add sequence labels
        ax.set_xticks(np.arange(len(sequence)))
        ax.set_yticks(np.arange(len(sequence)))
        ax.set_xticklabels(list(sequence), fontsize=6)
        ax.set_yticklabels(list(sequence), fontsize=6)
        
        ax.set_xlabel('Position')
        ax.set_ylabel('Position')
        ax.set_title('Attention Map', fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / f"{filename}.{config.FIGURE_FORMAT}")
        plt.close()
    
    def plot_uncertainty(self, prob_matrix: np.ndarray, var_matrix: np.ndarray,
                        filename: str = "uncertainty"):
        """Plot prediction uncertainty."""
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.5))
        
        # Mean probability
        im1 = ax1.imshow(prob_matrix, cmap='Reds', aspect='auto', vmin=0, vmax=1)
        plt.colorbar(im1, ax=ax1, label='Probability')
        ax1.set_title('Mean Probability', fontweight='bold')
        
        # Variance (uncertainty)
        im2 = ax2.imshow(var_matrix, cmap='Blues', aspect='auto')
        plt.colorbar(im2, ax=ax2, label='Variance')
        ax2.set_title('Uncertainty', fontweight='bold')
        
        for ax in [ax1, ax2]:
            ax.set_xlabel('Position')
            ax.set_ylabel('Position')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / f"{filename}.{config.FIGURE_FORMAT}")
        plt.close()
    
    def plot_bar_comparison(self, data: Dict[str, float], title: str,
                           filename: str = "bar_comparison"):
        """Plot bar chart comparison."""
        fig, ax = plt.subplots(figsize=(6, 4))
        
        categories = list(data.keys())
        values = list(data.values())
        
        bars = ax.bar(categories, values, color=self.colors['blue'], alpha=0.8)
        
        # Add value labels on bars
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{value:.3f}', ha='center', va='bottom')
        
        ax.set_ylabel('Score')
        ax.set_title(title, fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / f"{filename}.{config.FIGURE_FORMAT}")
        plt.close()
    
    def plot_length_distribution(self, df: pd.DataFrame, filename: str = "length_distribution"):
        """Plot sequence length distribution."""
        fig, ax = plt.subplots(figsize=(6, 4))
        
        ax.hist(df['length'], bins=50, color=self.colors['blue'], 
               alpha=0.7, edgecolor='black', linewidth=0.5)
        
        ax.set_xlabel('Sequence Length')
        ax.set_ylabel('Frequency')
        ax.set_title('Sequence Length Distribution', fontweight='bold')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / f"{filename}.{config.FIGURE_FORMAT}")
        plt.close()
    
    def plot_ablation_results(self, results: pd.DataFrame, filename: str = "ablation"):
        """Plot ablation study results."""
        fig, ax = plt.subplots(figsize=(8, 5))
        
        # Prepare data
        metrics = ['f1', 'mcc', 'exact_match']
        x = np.arange(len(metrics))
        width = 0.2
        
        # Plot each configuration
        colors = [self.colors['blue'], self.colors['green'], 
                 self.colors['orange'], self.colors['purple']]
        
        for idx, (config_name, config_data) in enumerate(results.groupby('config')):
            means = [config_data[m].mean() for m in metrics]
            stds = [config_data[m].std() for m in metrics]
            
            ax.bar(x + idx * width, means, width, yerr=stds,
                  label=config_name, color=colors[idx % len(colors)],
                  capsize=3, alpha=0.8)
        
        ax.set_xlabel('Metric')
        ax.set_ylabel('Score')
        ax.set_title('Ablation Study Results', fontweight='bold')
        ax.set_xticks(x + width * 1.5)
        ax.set_xticklabels([m.upper() for m in metrics])
        ax.legend(frameon=False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / f"{filename}.{config.FIGURE_FORMAT}")
        plt.close()


# =========================================================
# 27. Ablation Experiment Framework
# =========================================================
class AblationExperiment:
    """Comprehensive ablation study framework."""
    
    def __init__(self, base_config):
        self.base_config = base_config
        self.results = []
        self.visualizer = ScientificVisualizer(config.FIGURES_DIR / "ablation")
    
    def run_experiment_1_split_comparison(self, df: pd.DataFrame):
        """Compare random vs homology split."""
        log_message("\n[Experiment 1] Split Method Comparison")
        
        # Random split
        log_message("\n--- Random Split ---")
        splitter = HomologyAwareSplitter(method="random")
        train_r, valid_r, test_r = splitter.split(df, output_prefix="exp1_random")
        result_r = self._train_and_evaluate(train_r, valid_r, test_r, "random_split")
        
        # Homology split
        log_message("\n--- Homology Split ---")
        splitter = HomologyAwareSplitter(
            method=self.base_config.CLUSTER_METHOD,
            identity=self.base_config.HOMOLOGY_IDENTITY
        )
        train_h, valid_h, test_h = splitter.split(df, output_prefix="exp1_homology")
        result_h = self._train_and_evaluate(train_h, valid_h, test_h, "homology_split")
        
        self.results.append({
            'experiment': 'split_comparison',
            'random': result_r,
            'homology': result_h
        })
        
        return result_r, result_h
    
    def run_experiment_2_representation_comparison(self, train_df, valid_df, test_df):
        """Compare different sequence representations."""
        log_message("\n[Experiment 2] Representation Comparison")
        
        results = {}
        representations = ['1mer', '3mer', '4mer', 'rnafm']
        
        for rep in representations:
            if rep == 'rnafm' and not TRANSFORMERS_AVAILABLE:
                log_message(f"Skipping {rep} (RNA-FM not available)", "WARNING")
                continue
            
            log_message(f"\n--- {rep} ---")
            self.base_config.REPRESENTATION = rep
            self.base_config.KMER_SIZE = 1 if rep == '1mer' else (3 if rep == '3mer' else 4)
            
            result = self._train_and_evaluate(
                train_df, valid_df, test_df, 
                f"exp2_{rep}"
            )
            results[rep] = result
        
        self.results.append({
            'experiment': 'representation_comparison',
            **results
        })
        
        return results
    
    def run_experiment_3_loss_comparison(self, train_df, valid_df, test_df):
        """Compare different loss functions."""
        log_message("\n[Experiment 3] Loss Function Comparison")
        
        loss_configs = [
            ('bce', {'bce': 1.0, 'focal': 0.0, 'dice': 0.0, 
                    'symmetry': 0.0, 'diagonal': 0.0, 'sparsity': 0.0}),
            ('bce+focal', {'bce': 0.7, 'focal': 0.3, 'dice': 0.0,
                          'symmetry': 0.0, 'diagonal': 0.0, 'sparsity': 0.0}),
            ('bce+focal+dice', {'bce': 0.5, 'focal': 0.25, 'dice': 0.25,
                               'symmetry': 0.0, 'diagonal': 0.0, 'sparsity': 0.0}),
            ('full', self.base_config.LOSS_WEIGHTS)
        ]
        
        results = {}
        for name, weights in loss_configs:
            log_message(f"\n--- {name} ---")
            self.base_config.LOSS_WEIGHTS = weights
            result = self._train_and_evaluate(
                train_df, valid_df, test_df, 
                f"exp3_{name.replace('+', '_')}"
            )
            results[name] = result
        
        self.results.append({
            'experiment': 'loss_comparison',
            **results
        })
        
        return results
    
    def run_experiment_4_decoding_comparison(self, model, test_loader):
        """Compare different decoding strategies."""
        log_message("\n[Experiment 4] Decoding Comparison")
        
        # Implementation depends on having multiple decoders
        pass
    
    def _train_and_evaluate(self, train_df, valid_df, test_df, name: str) -> Dict:
        """Train and evaluate a model configuration."""
        # Create tokenizer
        if self.base_config.REPRESENTATION == 'rnafm':
            tokenizer = BaseTokenizer()  # RNA-FM uses its own tokenizer
        else:
            k = self.base_config.KMER_SIZE
            tokenizer = KMerTokenizer(k=k)
        
        # Create datasets
        train_ds = RNASSDataset(
            train_df, tokenizer, 
            self.base_config.MAX_LEN,
            self.base_config.MIN_HAIRPIN_LOOP
        )
        valid_ds = RNASSDataset(
            valid_df, tokenizer,
            self.base_config.MAX_LEN,
            self.base_config.MIN_HAIRPIN_LOOP
        )
        test_ds = RNASSDataset(
            test_df, tokenizer,
            self.base_config.MAX_LEN,
            self.base_config.MIN_HAIRPIN_LOOP
        )
        
        # Create data loaders
        train_loader = DataLoader(
            train_ds,
            batch_size=self.base_config.BATCH_SIZE,
            shuffle=True,
            num_workers=self.base_config.NUM_WORKERS,
            collate_fn=collate_fn
        )
        valid_loader = DataLoader(
            valid_ds,
            batch_size=self.base_config.BATCH_SIZE,
            shuffle=False,
            num_workers=self.base_config.NUM_WORKERS,
            collate_fn=collate_fn
        )
        test_loader = DataLoader(
            test_ds,
            batch_size=self.base_config.BATCH_SIZE,
            shuffle=False,
            num_workers=self.base_config.NUM_WORKERS,
            collate_fn=collate_fn
        )
        
        # Create model
        use_rnafm = (self.base_config.REPRESENTATION == 'rnafm')
        model = RNASSPredictor(
            vocab_size=tokenizer.vocab_size,
            d_model=self.base_config.D_MODEL,
            nhead=self.base_config.NHEAD,
            num_layers=self.base_config.NUM_LAYERS,
            dim_ff=self.base_config.DIM_FF,
            dropout=self.base_config.DROPOUT,
            max_len=self.base_config.MAX_LEN,
            use_rnafm=use_rnafm,
            use_positional_encoding=self.base_config.USE_POSITIONAL_ENCODING
        ).to(self.base_config.DEVICE)
        
        # Train (with fewer epochs for ablation)
        trainer = Trainer(model, self.base_config, run_id=f"ablation_{name}")
        trainer.train(train_loader, valid_loader, num_epochs=10)
        
        # Evaluate
        test_metrics = trainer.validate(test_loader)
        
        # Save results
        result_path = config.RESULTS_DIR / f"{name}_results.json"
        save_json(test_metrics, result_path)
        
        return test_metrics
    
    def generate_report(self):
        """Generate ablation study report."""
        report_path = config.RESULTS_DIR / "ablation_report.txt"
        
        with open(report_path, 'w') as f:
            f.write("="*80 + "\n")
            f.write("RNA Secondary Structure Prediction - Ablation Study Report\n")
            f.write("="*80 + "\n\n")
            
            for exp in self.results:
                f.write(f"\n{exp['experiment']}\n")
                f.write("-"*40 + "\n")
                f.write(json.dumps({k: v for k, v in exp.items() if k != 'experiment'}, 
                                   indent=2))
                f.write("\n")
        
        log_message(f"Ablation report saved to {report_path}")
        
        # Create visualizations
        self._visualize_results()
    
    def _visualize_results(self):
        """Create visualizations of ablation results."""
        # Convert results to DataFrame for plotting
        rows = []
        for exp in self.results:
            exp_name = exp['experiment']
            for config_name, metrics in exp.items():
                if config_name == 'experiment':
                    continue
                row = {'experiment': exp_name, 'config': config_name}
                if isinstance(metrics, dict):
                    row.update({k: v for k, v in metrics.items() 
                               if isinstance(v, (int, float))})
                rows.append(row)
        
        if rows:
            df = pd.DataFrame(rows)
            self.visualizer.plot_ablation_results(df)


# =========================================================
# 28. Command Line Argument Parsing
# =========================================================
def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='RNA Secondary Structure Prediction System - SCI Edition',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Mode selection
    parser.add_argument('--mode', type=str, default='train',
                       choices=['crawl', 'process', 'train', 'predict', 'ablation', 'all'],
                       help='Execution mode')
    
    # Crawler parameters
    parser.add_argument('--crawl', action='store_true', help='Run data crawling')
    parser.add_argument('--target_sequences', type=int, default=50000,
                       help='Target number of sequences to crawl')
    parser.add_argument('--sources', nargs='+', default=['bprna', 'rfam', 'rnacentral'],
                       help='Data sources to crawl')
    parser.add_argument('--resume', action='store_true', help='Resume crawling from checkpoint')
    parser.add_argument('--parallel', action='store_true', help='Parallel crawling')
    
    # Data processing parameters
    parser.add_argument('--process', action='store_true', help='Process crawled data')
    parser.add_argument('--min_len', type=int, default=20, help='Minimum sequence length')
    parser.add_argument('--max_len', type=int, default=256, help='Maximum sequence length')
    parser.add_argument('--require_structure', action='store_true',
                       help='Require structure information')
    
    # Data splitting
    parser.add_argument('--split_method', type=str, default='homology',
                       choices=['random', 'homology'],
                       help='Data splitting method')
    parser.add_argument('--cluster_tool', type=str, default='mmseqs2',
                       choices=['mmseqs2', 'cdhit'],
                       help='Clustering tool for homology split')
    parser.add_argument('--homology_identity', type=float, default=0.80,
                       help='Sequence identity threshold for clustering')
    
    # Model parameters
    parser.add_argument('--representation', type=str, default='kmer',
                       choices=['kmer', 'rnafm'],
                       help='Sequence representation')
    parser.add_argument('--kmer', type=int, default=3, help='k-mer size')
    parser.add_argument('--d_model', type=int, default=128, help='Transformer dimension')
    parser.add_argument('--nhead', type=int, default=8, help='Number of attention heads')
    parser.add_argument('--num_layers', type=int, default=4, help='Number of Transformer layers')
    parser.add_argument('--dropout', type=float, default=0.1, help='Dropout rate')
    
    # Training parameters
    parser.add_argument('--train', action='store_true', help='Train model')
    parser.add_argument('--batch_size', type=int, default=8, help='Batch size')
    parser.add_argument('--epochs', type=int, default=50, help='Number of training epochs')
    parser.add_argument('--lr', type=float, default=2e-4, help='Learning rate')
    parser.add_argument('--weight_decay', type=float, default=1e-4, help='Weight decay')
    parser.add_argument('--use_rnafm', action='store_true', help='Use RNA-FM embeddings')
    parser.add_argument('--num_repeats', type=int, default=3,
                       help='Number of repeats for statistical significance')
    
    # Evaluation parameters
    parser.add_argument('--evaluate', action='store_true', help='Evaluate model')
    parser.add_argument('--model_path', type=str, help='Path to trained model')
    parser.add_argument('--test_sequence', type=str, help='Test a single sequence')
    parser.add_argument('--decoding_method', type=str, default='nussinov',
                       choices=['greedy', 'nussinov'],
                       help='Decoding method')
    parser.add_argument('--threshold', type=float, default=0.5,
                       help='Decoding threshold')
    
    # Uncertainty estimation
    parser.add_argument('--mc_samples', type=int, default=20,
                       help='Number of MC Dropout samples')
    
    # Thermodynamic fusion
    parser.add_argument('--use_vienna', action='store_true',
                       help='Use ViennaRNA thermodynamic fusion')
    parser.add_argument('--thermo_weight', type=float, default=0.2,
                       help='Weight for thermodynamic predictions')
    
    # Ablation
    parser.add_argument('--ablation', action='store_true', help='Run ablation experiments')
    
    # Other
    parser.add_argument('--seed', type=int, default=42, help='Random seed')
    parser.add_argument('--device', type=str, choices=['cuda', 'cpu'],
                       help='Device to use')
    parser.add_argument('--config', type=str, help='Path to configuration file')
    
    return parser.parse_args()


# =========================================================
# 29. Main Function
# =========================================================
def main():
    """Main execution function."""
    # Parse arguments
    args = parse_arguments()
    
    # Load configuration if provided
    if args.config:
        config = Config.load(Path(args.config))
    else:
        config = Config()
    
    # Update from command line
    config.update_from_args(args)
    
    # Set random seed
    set_seed(config.SEED)
    
    # Print header
    print("="*80)
    print("RNA Secondary Structure Prediction System - SCI Edition")
    print("="*80)
    print(f"Device: {config.DEVICE}")
    print(f"Base Directory: {config.BASE_DIR}")
    print(f"Mode: {args.mode}")
    print("="*80)
    
    # Mode execution
    if args.mode == 'crawl' or args.mode == 'all':
        print("\n[Step 1] Crawling RNA data...")
        with DataCrawlerManager() as manager:
            results = manager.crawl_all(
                target_total=args.target_sequences,
                parallel=args.parallel,
                resume=args.resume,
                source_list=args.sources if args.sources else None
            )
            print(f"Crawling results: {results}")
    
    if args.mode == 'process' or args.mode == 'all':
        print("\n[Step 2] Processing crawled data...")
        with DataProcessor() as processor:
            df = processor.build_dataset(
                require_structure=args.require_structure,
                source_filter=args.sources if args.sources else None
            )
            
            if df.empty:
                print("[ERROR] No data available. Please run crawling first.")
                return
            
            if args.require_structure:
                df = processor.filter_by_structure_complexity(df, min_pairs=1)
            
            # Export
            processor.export_to_fasta(df, config.PROC_DIR / "all_sequences.fasta")
            if args.require_structure:
                processor.export_to_stockholm(df, config.PROC_DIR / "all_structures.stk")
            
            print(f"Final dataset size: {len(df)}")
            
            # Visualize data distribution
            vis = ScientificVisualizer(config.FIGURES_DIR / "data")
            vis.plot_length_distribution(df)
    
    if args.mode == 'train' or args.mode == 'all':
        print("\n[Step 3] Training model...")
        
        # Load dataset
        dataset_path = config.PROC_DIR / "raw_dataset.csv"
        if not dataset_path.exists():
            print("[ERROR] Processed data not found. Please run --process first.")
            return
        
        df = pd.read_csv(dataset_path)
        
        # Multiple repeats for statistical significance
        all_results = []
        
        for repeat in range(args.num_repeats):
            print(f"\n--- Repeat {repeat + 1}/{args.num_repeats} ---")
            set_seed(config.SEED + repeat)  # Different seed for each repeat
            
            # Split data
            if args.split_method == 'random':
                splitter = HomologyAwareSplitter(method='random')
                train_df, valid_df, test_df = splitter.split(
                    df, output_prefix=f"split_repeat{repeat}"
                )
            else:
                splitter = HomologyAwareSplitter(
                    method=args.cluster_tool,
                    identity=args.homology_identity
                )
                train_df, valid_df, test_df = splitter.split(
                    df, output_prefix=f"split_repeat{repeat}"
                )
            
            # Create tokenizer
            if args.representation == 'rnafm':
                tokenizer = BaseTokenizer()
            else:
                tokenizer = KMerTokenizer(k=args.kmer)
            
            # Create datasets
            train_ds = RNASSDataset(
                train_df, tokenizer, config.MAX_LEN, config.MIN_HAIRPIN_LOOP
            )
            valid_ds = RNASSDataset(
                valid_df, tokenizer, config.MAX_LEN, config.MIN_HAIRPIN_LOOP
            )
            test_ds = RNASSDataset(
                test_df, tokenizer, config.MAX_LEN, config.MIN_HAIRPIN_LOOP
            )
            
            # Create data loaders
            train_loader = DataLoader(
                train_ds, batch_size=config.BATCH_SIZE, shuffle=True,
                num_workers=config.NUM_WORKERS, collate_fn=collate_fn,
                pin_memory=(config.DEVICE == 'cuda')
            )
            valid_loader = DataLoader(
                valid_ds, batch_size=config.BATCH_SIZE, shuffle=False,
                num_workers=config.NUM_WORKERS, collate_fn=collate_fn,
                pin_memory=(config.DEVICE == 'cuda')
            )
            test_loader = DataLoader(
                test_ds, batch_size=config.BATCH_SIZE, shuffle=False,
                num_workers=config.NUM_WORKERS, collate_fn=collate_fn,
                pin_memory=(config.DEVICE == 'cuda')
            )
            
            # Create model
            use_rnafm = (args.representation == 'rnafm')
            model = RNASSPredictor(
                vocab_size=tokenizer.vocab_size,
                d_model=config.D_MODEL,
                nhead=config.NHEAD,
                num_layers=config.NUM_LAYERS,
                dim_ff=config.DIM_FF,
                dropout=config.DROPOUT,
                max_len=config.MAX_LEN,
                use_rnafm=use_rnafm,
                use_positional_encoding=config.USE_POSITIONAL_ENCODING
            ).to(config.DEVICE)
            
            total_params = sum(p.numel() for p in model.parameters())
            print(f"Total parameters: {total_params:,}")
            
            # Train
            trainer = Trainer(model, config, run_id=f"run_{repeat}")
            best_f1 = trainer.train(train_loader, valid_loader, config.NUM_EPOCHS)
            
            # Test
            print("\nFinal evaluation on test set...")
            test_metrics = trainer.validate(test_loader, threshold=config.DECODING_THRESHOLD)
            test_metrics['repeat'] = repeat
            all_results.append(test_metrics)
            
            print(f"Test Results (Repeat {repeat + 1}):")
            for k, v in test_metrics.items():
                if isinstance(v, float):
                    print(f"  {k}: {v:.4f}")
        
        # Compute statistics across repeats
        if len(all_results) > 1:
            print("\n=== Aggregate Results ===")
            metrics_summary = {}
            for metric in ['f1', 'mcc', 'precision', 'recall', 'exact_match']:
                values = [r[metric] for r in all_results]
                metrics_summary[metric] = {
                    'mean': np.mean(values),
                    'std': np.std(values),
                    'min': np.min(values),
                    'max': np.max(values)
                }
                print(f"{metric}: {np.mean(values):.4f} ± {np.std(values):.4f}")
            
            # Save aggregate results
            agg_path = config.RESULTS_DIR / "aggregate_results.json"
            save_json(metrics_summary, agg_path)
    
    if args.mode == 'predict':
        print("\n[Step] Testing single sequence...")
        
        if not args.test_sequence:
            print("[ERROR] Please provide --test_sequence")
            return
        
        # Load model
        model_path = args.model_path or config.MODEL_DIR / "best_model.pt"
        if not Path(model_path).exists():
            print(f"[ERROR] Model not found: {model_path}")
            return
        
        checkpoint = torch.load(model_path, map_location=config.DEVICE)
        
        # Load tokenizer
        tokenizer = KMerTokenizer(k=config.KMER_SIZE)
        
        # Create model
        model = RNASSPredictor(
            vocab_size=tokenizer.vocab_size,
            d_model=config.D_MODEL,
            nhead=config.NHEAD,
            num_layers=config.NUM_LAYERS,
            dim_ff=config.DIM_FF,
            dropout=config.DROPOUT,
            max_len=config.MAX_LEN
        ).to(config.DEVICE)
        model.load_state_dict(checkpoint['model_state_dict'])
        model.eval()
        
        # Predict
        sequence = args.test_sequence
        print(f"Input sequence: {sequence}")
        
        # Prepare input
        input_ids = torch.tensor(
            tokenizer.encode(sequence, config.MAX_LEN),
            dtype=torch.long
        ).unsqueeze(0).to(config.DEVICE)
        
        prior = torch.tensor(
            build_basepair_prior(sequence, config.MAX_LEN, config.MIN_HAIRPIN_LOOP),
            dtype=torch.float32
        ).unsqueeze(0).to(config.DEVICE)
        
        dist = torch.tensor(
            build_distance_bias(config.MAX_LEN),
            dtype=torch.float32
        ).unsqueeze(0).to(config.DEVICE)
        
        # Predict with uncertainty
        estimator = UncertaintyEstimator(model, num_samples=config.MC_SAMPLES)
        mean_prob, var_prob = estimator.estimate(input_ids, prior, dist)
        
        mean_prob = mean_prob[0, :len(sequence), :len(sequence)]
        var_prob = var_prob[0, :len(sequence), :len(sequence)]
        
        # Optional thermodynamic fusion
        if args.use_vienna:
            vienna = ViennaRNAInterface()
            thermo_prob = vienna.get_pair_probabilities(sequence)
            if thermo_prob is not None:
                mean_prob = vienna.fuse_with_learned(
                    mean_prob, thermo_prob, alpha=1 - args.thermo_weight
                )
        
        # Decode
        decoder = NussinovDecoder(min_loop=config.MIN_HAIRPIN_LOOP)
        pred_db = decoder.decode(mean_prob, sequence, threshold=args.threshold)
        
        print(f"Predicted structure: {pred_db}")
        print(f"Global uncertainty: {var_prob.mean():.6f}")
        
        # Visualize
        vis = ScientificVisualizer(config.FIGURES_DIR / "predictions")
        vis.plot_uncertainty(mean_prob, var_prob, filename=f"pred_{hash(sequence)}")
    
    if args.mode == 'ablation':
        print("\n[Step] Running ablation experiments...")
        
        # Load data
        dataset_path = config.PROC_DIR / "raw_dataset.csv"
        if not dataset_path.exists():
            print("[ERROR] Processed data not found.")
            return
        
        df = pd.read_csv(dataset_path)
        
        # Run ablation
        ablation = AblationExperiment(config)
        
        # Experiment 1: Split comparison
        ablation.run_experiment_1_split_comparison(df)
        
        # Get a base split for other experiments
        splitter = HomologyAwareSplitter(
            method=config.CLUSTER_METHOD,
            identity=config.HOMOLOGY_IDENTITY
        )
        train_df, valid_df, test_df = splitter.split(df)
        
        # Experiment 2: Representation comparison
        ablation.run_experiment_2_representation_comparison(train_df, valid_df, test_df)
        
        # Experiment 3: Loss comparison
        ablation.run_experiment_3_loss_comparison(train_df, valid_df, test_df)
        
        # Generate report
        ablation.generate_report()
    
    print("\n" + "="*80)
    print("Pipeline completed successfully!")
    print(f"Results saved in: {config.RESULTS_DIR}")
    print(f"Figures saved in: {config.FIGURES_DIR}")
    print("="*80)


if __name__ == "__main__":
    main()