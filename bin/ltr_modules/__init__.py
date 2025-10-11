#!/usr/bin/env python3
"""
LTR Modules Package

Modular components for LTR retrotransposon analysis with better
code organization, maintainability, and reusability.

Modules:
- LTR_QualityValidator: Quality validation and scoring
- LTR_BoundaryRefiner: Boundary detection and extension
- LTR_ChimeraDetector: Chimera detection and splitting
- LTR_ConsensusBuilder: Consensus sequence building
- shared_utils: Shared utilities and constants
- genome_utils: Genome access utilities
"""

from .LTR_QualityValidator import LTRQualityValidator
from .LTR_BoundaryRefiner import LTRBoundaryRefiner
from .LTR_ChimeraDetector import LTRChimeraDetector
from .LTR_ConsensusBuilder import LTRConsensusBuilder
from .genome_utils import GenomeAccess, GenomeSeqkitAccess, GenomeSimpleAccess
from .shared_utils import (
    LTRConstants,
    calculate_similarity,
    reverse_complement,
    gc_content,
    check_low_complexity,
    find_tandem_repeats,
    validate_sequence,
    compute_pairwise_similarity,
    compute_avg_similarity
)

__all__ = [
    # Main classes
    'LTRQualityValidator',
    'LTRBoundaryRefiner',
    'LTRChimeraDetector',
    'LTRConsensusBuilder',
    # Genome access
    'GenomeAccess',
    'GenomeSeqkitAccess',
    'GenomeSimpleAccess',
    # Shared utilities
    'LTRConstants',
    'calculate_similarity',
    'reverse_complement',
    'gc_content',
    'check_low_complexity',
    'find_tandem_repeats',
    'validate_sequence',
    'compute_pairwise_similarity',
    'compute_avg_similarity',
]

__version__ = '1.0.0'
__author__ = 'Pan_TE Development Team'
__date__ = '2025-10-07'
