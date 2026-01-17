from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Settings:
    # Runtime inputs
    tasks: tuple[str, ...]
    ids_path: Path
    out_dir: Path
    batch_size: int
    output_format: str

    comments: CommentsTemplates | None = None
    cross_references: CrossReferencesTemplates | None = None
    ec: ECTemplates | None = None
    features: FeaturesTemplates | None = None

    @property
    def _out_suffix(self):
        if self.output_format in ("jsonl", "csv", "tsv"):
            return f".{self.output_format}"
        raise ValueError(f"Unsupported output_format: {self.output_format}")

    # comments outputs
    @property
    def ca_out(self):
        return self.out_dir / "catalytic_activity" / f"catalytic_activity{self._out_suffix}"

    @property
    def cf_out(self):
        return self.out_dir / "cofactor" / f"cofactor{self._out_suffix}"

    @property
    def sl_out(self):
        return self.out_dir / "subcellular_location" / f"subcellular_location{self._out_suffix}"

    @property
    def pw_out(self):
        return self.out_dir / "pathway" / f"pathway{self._out_suffix}"

    # cross_references outputs
    @property
    def go_bp_out(self):
        return self.out_dir / "gene_ontology" / f"go_bp{self._out_suffix}"

    @property
    def go_mf_out(self):
        return self.out_dir / "gene_ontology" / f"go_mf{self._out_suffix}"

    @property
    def go_cc_out(self):
        return self.out_dir / "gene_ontology" / f"go_cc{self._out_suffix}"

    # ec outputs
    @property
    def ec_out(self):
        return self.out_dir / "ec" / f"ec{self._out_suffix}"

    # features outputs
    @property
    def rg_out(self):
        return self.out_dir / "regions" / f"regions{self._out_suffix}"

    @property
    def st_out(self):
        return self.out_dir / "sites" / f"sites{self._out_suffix}"

    def __post_init__(self):
        out_files = []
        tasks_set = set(self.tasks)

        if "comments" in tasks_set:
            out_files.extend([self.ca_out, self.cf_out, self.sl_out, self.pw_out])

        if "cross_references" in tasks_set:
            out_files.extend([self.go_bp_out, self.go_mf_out, self.go_cc_out])

        if "ec" in tasks_set:
            out_files.extend([self.ec_out])

        if "features" in tasks_set:
            out_files.extend([self.rg_out, self.st_out])

        if not out_files:
            raise ValueError("No tasks selected: Settings.tasks is empty")

        for p in out_files:
            p.parent.mkdir(parents=True, exist_ok=True)

    def validate(self):
        if not self.ids_path.exists():
            raise FileNotFoundError(f"IDs file not found: {self.ids_path}")


@dataclass(frozen=True)
class CommentsTemplates:
    ca_a: Path | None = None
    cf_a: Path | None = None
    cf_b: Path | None = None
    sl_a: Path | None = None
    sl_b: Path | None = None
    sl_c: Path | None = None
    sl_d: Path | None = None
    pw_a: Path | None = None
    pw_b: Path | None = None
    pw_c: Path | None = None
    pw_d: Path | None = None


@dataclass(frozen=True)
class CrossReferencesTemplates:
    go_bp: Path | None = None
    go_mf: Path | None = None
    go_cc: Path | None = None


@dataclass(frozen=True)
class ECTemplates:
    ec_a: Path | None = None
    ec_b: Path | None = None
    ec_c: Path | None = None
    ec_d: Path | None = None
    ec_e: Path | None = None
    ec_f: Path | None = None
    ec_g: Path | None = None
    ec_h: Path | None = None
    ec_i: Path | None = None
    ec_j: Path | None = None
    ec_k: Path | None = None


@dataclass(frozen=True)
class FeaturesTemplates:
    rg_a: Path | None = None
    rg_b: Path | None = None
    st_a: Path | None = None
    st_b: Path | None = None
    st_c: Path | None = None
    st_d: Path | None = None