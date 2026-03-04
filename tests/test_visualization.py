"""Tests for uropa/visualization.py"""

import os
import logging
import tempfile

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pytest

from uropa.visualization import (
    load_uropa_output,
    reconstruct_signed_distance,
    plot_distance_kde,
    plot_distance_hist,
    generate_summary,
)


# ── Fixtures ──────────────────────────────────────────────────────────────────


@pytest.fixture
def sample_df():
    """Minimal DataFrame mimicking UROPA output structure."""
    return pd.DataFrame({
        'peak_chr': ['chr1'] * 6,
        'peak_start': [100, 200, 300, 100, 200, 300],
        'peak_end': [200, 300, 400, 200, 300, 400],
        'feature': ['gene', 'gene', 'transcript', 'gene', 'gene', 'transcript'],
        'feat_start': [50, 150, 250, 50, 150, 250],
        'feat_end': [250, 350, 450, 250, 350, 450],
        'feat_strand': ['+', '+', '+', '-', '-', '-'],
        'feat_anchor': ['start', 'end', 'center', 'start', 'end', 'center'],
        'distance': [100, 150, 0, 80, 200, 50],
        'name': ['q1', 'q1', 'q1', 'q2', 'q2', 'q2'],
    })


@pytest.fixture
def reference_prefix():
    """Path prefix for reference output files."""
    return os.path.join(os.path.dirname(__file__), '..', 'reference_output', 'genomic_regions')


@pytest.fixture
def logger():
    """Simple logger for testing."""
    log = logging.getLogger('test_visualization')
    log.setLevel(logging.DEBUG)
    if not log.handlers:
        log.addHandler(logging.StreamHandler())
    return log


# ── load_uropa_output ─────────────────────────────────────────────────────────


class TestLoadUropaOutput:

    def test_loads_valid_tsv(self, tmp_path):
        tsv = tmp_path / "test.txt"
        tsv.write_text(
            "peak_start\tpeak_end\tfeat_start\tfeat_end\tdistance\n"
            "100\t200\t50\t250\t100\n"
            "300\t400\t250\t450\t0\n"
        )
        df = load_uropa_output(str(tsv))
        assert len(df) == 2
        assert df['distance'].dtype == float

    def test_drops_na_rows(self, tmp_path):
        tsv = tmp_path / "test.txt"
        tsv.write_text(
            "peak_start\tpeak_end\tfeat_start\tfeat_end\tdistance\n"
            "100\t200\t50\t250\t100\n"
            "300\t400\t250\t450\tNA\n"
        )
        df = load_uropa_output(str(tsv))
        assert len(df) == 1

    def test_distance_is_float(self, tmp_path):
        tsv = tmp_path / "test.txt"
        tsv.write_text(
            "peak_start\tpeak_end\tfeat_start\tfeat_end\tdistance\n"
            "100\t200\t50\t250\t42\n"
        )
        df = load_uropa_output(str(tsv))
        assert df['distance'].dtype == float


# ── reconstruct_signed_distance ──────────────────────────────────────────────


class TestReconstructSignedDistance:

    def _make_row(self, feat_strand, feat_anchor):
        """Helper: single-row df with peak_center=150, feat 100-200."""
        return pd.DataFrame({
            'peak_start': [100.0],
            'peak_end': [200.0],
            'feat_start': [100.0],
            'feat_end': [200.0],
            'feat_strand': [feat_strand],
            'feat_anchor': [feat_anchor],
        })

    def test_plus_strand_start_anchor(self):
        # anchor=feat_start=100, peak_center=150, raw=50, signed=+50
        df = reconstruct_signed_distance(self._make_row('+', 'start'))
        assert df['signed_distance'].iloc[0] == pytest.approx(50.0)

    def test_minus_strand_start_anchor(self):
        # anchor=feat_end=200 (swapped), peak_center=150, raw=-50, signed=+50 (flipped)
        df = reconstruct_signed_distance(self._make_row('-', 'start'))
        assert df['signed_distance'].iloc[0] == pytest.approx(50.0)

    def test_plus_strand_end_anchor(self):
        # anchor=feat_end=200, peak_center=150, raw=-50, signed=-50
        df = reconstruct_signed_distance(self._make_row('+', 'end'))
        assert df['signed_distance'].iloc[0] == pytest.approx(-50.0)

    def test_minus_strand_end_anchor(self):
        # anchor=feat_start=100 (swapped), peak_center=150, raw=50, signed=-50 (flipped)
        df = reconstruct_signed_distance(self._make_row('-', 'end'))
        assert df['signed_distance'].iloc[0] == pytest.approx(-50.0)

    def test_center_anchor_plus(self):
        # anchor=midpoint=150, peak_center=150, raw=0, signed=0
        df = reconstruct_signed_distance(self._make_row('+', 'center'))
        assert df['signed_distance'].iloc[0] == pytest.approx(0.0)

    def test_center_anchor_minus(self):
        # anchor=midpoint=150, peak_center=150, raw=0, signed=0
        df = reconstruct_signed_distance(self._make_row('-', 'center'))
        assert df['signed_distance'].iloc[0] == pytest.approx(0.0)


# ── Plot functions ────────────────────────────────────────────────────────────


class TestPlotDistanceKde:

    def test_creates_plot_elements(self, sample_df):
        fig, ax = plt.subplots()
        plot_distance_kde(sample_df, 'feature', 'distance', False, ax)
        # Should have lines and collections (fill_between) on axes
        assert len(ax.lines) > 0 or len(ax.collections) > 0
        plt.close(fig)

    def test_no_crash_with_single_point(self):
        """Single data point per group should be skipped, not crash."""
        df = pd.DataFrame({
            'feature': ['gene'],
            'distance': [100.0],
        })
        fig, ax = plt.subplots()
        plot_distance_kde(df, 'feature', 'distance', False, ax)
        plt.close(fig)


class TestPlotDistanceHist:

    def test_creates_patches(self, sample_df):
        fig, ax = plt.subplots()
        plot_distance_hist(sample_df, 'feature', 'distance', False, ax)
        assert len(ax.patches) > 0
        plt.close(fig)


# ── generate_summary (integration) ───────────────────────────────────────────


class TestGenerateSummary:

    def test_finalhits_kde_directional(self, reference_prefix, logger, tmp_path):
        """Default params: finalhits, kde, directional, feature facet."""
        prefix = str(tmp_path / 'genomic_regions')
        # Symlink reference files into tmp_path
        os.symlink(os.path.abspath(reference_prefix + '_finalhits.txt'),
                    prefix + '_finalhits.txt')
        generate_summary(prefix, 'finalhits', 'kde', 'feature',
                         'directional', False, logger)
        pdf = prefix + '_summary.pdf'
        assert os.path.exists(pdf)
        assert os.path.getsize(pdf) > 0

    def test_finalhits_hist_absolute(self, reference_prefix, logger, tmp_path):
        prefix = str(tmp_path / 'genomic_regions')
        os.symlink(os.path.abspath(reference_prefix + '_finalhits.txt'),
                    prefix + '_finalhits.txt')
        generate_summary(prefix, 'finalhits', 'hist', 'feature',
                         'absolute', False, logger)
        pdf = prefix + '_summary.pdf'
        assert os.path.exists(pdf)
        assert os.path.getsize(pdf) > 0

    def test_allhits_kde(self, reference_prefix, logger, tmp_path):
        prefix = str(tmp_path / 'genomic_regions')
        os.symlink(os.path.abspath(reference_prefix + '_allhits.txt'),
                    prefix + '_allhits.txt')
        generate_summary(prefix, 'allhits', 'kde', 'feature',
                         'directional', False, logger)
        pdf = prefix + '_summary.pdf'
        assert os.path.exists(pdf)
        assert os.path.getsize(pdf) > 0

    def test_log_scale(self, reference_prefix, logger, tmp_path):
        prefix = str(tmp_path / 'genomic_regions')
        os.symlink(os.path.abspath(reference_prefix + '_finalhits.txt'),
                    prefix + '_finalhits.txt')
        generate_summary(prefix, 'finalhits', 'kde', 'feature',
                         'directional', True, logger)
        pdf = prefix + '_summary.pdf'
        assert os.path.exists(pdf)


# ── Edge cases ────────────────────────────────────────────────────────────────


class TestEdgeCases:

    def test_missing_file(self, logger, tmp_path):
        """Missing input file logs warning, no crash."""
        prefix = str(tmp_path / 'nonexistent')
        generate_summary(prefix, 'finalhits', 'kde', 'feature',
                         'directional', False, logger)
        # No exception, no PDF created
        assert not os.path.exists(prefix + '_summary.pdf')

    def test_all_na_distances(self, logger, tmp_path):
        """All NA distances: returns without error, no PDF."""
        prefix = str(tmp_path / 'empty')
        tsv = prefix + '_finalhits.txt'
        with open(tsv, 'w') as f:
            f.write("peak_start\tpeak_end\tfeat_start\tfeat_end\tdistance\tfeature\t"
                    "feat_strand\tfeat_anchor\tname\n")
            f.write("100\t200\t50\t250\tNA\tgene\t+\tstart\tq1\n")
        generate_summary(prefix, 'finalhits', 'kde', 'feature',
                         'directional', False, logger)
        assert not os.path.exists(prefix + '_summary.pdf')
