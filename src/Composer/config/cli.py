import argparse
from pathlib import Path
from .settings import Settings, CommentsTemplates, CrossReferencesTemplates, ECTemplates, FeaturesTemplates


def build_parser():
    p = argparse.ArgumentParser()

    # Runtime inputs
    p.add_argument(
        "--tasks",
        required=True,
        nargs="+",
        choices=("comments", "features", "ec", "cross_references"),
        help='Which task(s) to run (space-separated). Choices: comments | features | ec | cross_references',
    )
    p.add_argument("--ids", required=True, help="Protein ID list file (txt)")
    p.add_argument("--out_dir", required=True, help="Output directory")
    p.add_argument("--batch_size", type=int, default=100, help="Write buffer size")
    p.add_argument("--output_format", default="jsonl", help="Output format: jsonl | csv | tsv")

    # comments
    p.add_argument("--ca_a", default=None, help="Catalytic activity template file A (json)")
    p.add_argument("--cf_a", default=None, help="Cofactor template file A (json)")
    p.add_argument("--cf_b", default=None, help="Cofactor template file B (json)")
    p.add_argument("--sl_a", default=None, help="Subcellular location template file A (json)")
    p.add_argument("--sl_b", default=None, help="Subcellular location template file B (json)")
    p.add_argument("--sl_c", default=None, help="Subcellular location template file C (json)")
    p.add_argument("--sl_d", default=None, help="Subcellular location template file D (json)")
    p.add_argument("--pw_a", default=None, help="Pathway template file A (json)")
    p.add_argument("--pw_b", default=None, help="Pathway template file B (json)")
    p.add_argument("--pw_c", default=None, help="Pathway template file C (json)")
    p.add_argument("--pw_d", default=None, help="Pathway template file D (json)")

    # cross_references
    p.add_argument("--go_bp", default=None, help="GO Biological Process templates (json)")
    p.add_argument("--go_mf", default=None, help="GO Molecular Function templates (json)")
    p.add_argument("--go_cc", default=None, help="GO Cellular Component templates (json)")

    # ec
    p.add_argument("--ec_a", default=None, help="EC templates A (json)")
    p.add_argument("--ec_b", default=None, help="EC templates B (json)")
    p.add_argument("--ec_c", default=None, help="EC templates C (json)")
    p.add_argument("--ec_d", default=None, help="EC templates D (json)")
    p.add_argument("--ec_e", default=None, help="EC templates E (json)")
    p.add_argument("--ec_f", default=None, help="EC templates F (json)")
    p.add_argument("--ec_g", default=None, help="EC templates G (json)")
    p.add_argument("--ec_h", default=None, help="EC templates H (json)")
    p.add_argument("--ec_i", default=None, help="EC templates I (json)")
    p.add_argument("--ec_j", default=None, help="EC templates J (json)")
    p.add_argument("--ec_k", default=None, help="EC templates K (json)")

    # features
    p.add_argument("--rg_a", default=None, help="Regions templates A (json)")
    p.add_argument("--rg_b", default=None, help="Regions templates B (json)")
    p.add_argument("--st_a", default=None, help="Sites templates A (json)")
    p.add_argument("--st_b", default=None, help="Sites templates B (json)")
    p.add_argument("--st_c", default=None, help="Sites templates C (json)")
    p.add_argument("--st_d", default=None, help="Sites templates D (json)")

    return p


def _to_path(v):
    return Path(v) if v is not None else None


def _validate_required(task_name, required):
    missing = [k for k, v in required.items() if v is None]
    if missing:
        raise ValueError(f"task={task_name} missing template arg(s): {', '.join(missing)}")

    not_found = [k for k, p in required.items() if p is not None and not p.exists()]
    if not_found:
        raise FileNotFoundError(
            f"task={task_name} template file(s) not found: "
            + ", ".join(f"{k}={required[k]}" for k in not_found)
        )


def parse_settings():
    args = build_parser().parse_args()

    tasks = tuple(args.tasks)
    tasks_set = set(tasks)

    comments = None
    cross_references = None
    ec = None
    features = None

    if "comments" in tasks_set:
        comments = CommentsTemplates(
            ca_a=_to_path(args.ca_a),
            cf_a=_to_path(args.cf_a),
            cf_b=_to_path(args.cf_b),
            sl_a=_to_path(args.sl_a),
            sl_b=_to_path(args.sl_b),
            sl_c=_to_path(args.sl_c),
            sl_d=_to_path(args.sl_d),
            pw_a=_to_path(args.pw_a),
            pw_b=_to_path(args.pw_b),
            pw_c=_to_path(args.pw_c),
            pw_d=_to_path(args.pw_d),
        )

        t = comments
        required = {
            "ca_a": t.ca_a,
            "cf_a": t.cf_a, "cf_b": t.cf_b,
            "sl_a": t.sl_a, "sl_b": t.sl_b, "sl_c": t.sl_c, "sl_d": t.sl_d,
            "pw_a": t.pw_a, "pw_b": t.pw_b, "pw_c": t.pw_c, "pw_d": t.pw_d,
        }
        _validate_required("comments", required)

    if "cross_references" in tasks_set:
        cross_references = CrossReferencesTemplates(
            go_bp=_to_path(args.go_bp),
            go_mf=_to_path(args.go_mf),
            go_cc=_to_path(args.go_cc),
        )

        t = cross_references
        required = {
            "go_bp": t.go_bp, "go_mf": t.go_mf, "go_cc": t.go_cc,
        }
        _validate_required("cross_references", required)

    if "ec" in tasks_set:
        ec = ECTemplates(
            ec_a=_to_path(args.ec_a),
            ec_b=_to_path(args.ec_b),
            ec_c=_to_path(args.ec_c),
            ec_d=_to_path(args.ec_d),
            ec_e=_to_path(args.ec_e),
            ec_f=_to_path(args.ec_f),
            ec_g=_to_path(args.ec_g),
            ec_h=_to_path(args.ec_h),
            ec_i=_to_path(args.ec_i),
            ec_j=_to_path(args.ec_j),
            ec_k=_to_path(args.ec_k),
        )

        t = ec
        required = {
            "ec_a": t.ec_a, "ec_b": t.ec_b, "ec_c": t.ec_c, "ec_d": t.ec_d,
            "ec_e": t.ec_e, "ec_f": t.ec_f, "ec_g": t.ec_g, "ec_h": t.ec_h,
            "ec_i": t.ec_i, "ec_j": t.ec_j, "ec_k": t.ec_k,
        }
        _validate_required("ec", required)

    if "features" in tasks_set:
        features = FeaturesTemplates(
            rg_a=_to_path(args.rg_a),
            rg_b=_to_path(args.rg_b),
            st_a=_to_path(args.st_a),
            st_b=_to_path(args.st_b),
            st_c=_to_path(args.st_c),
            st_d=_to_path(args.st_d),
        )

        t = features
        required = {
            "rg_a": t.rg_a, "rg_b": t.rg_b,
            "st_a": t.st_a, "st_b": t.st_b, "st_c": t.st_c, "st_d": t.st_d,
        }
        _validate_required("features", required)

    s = Settings(
        tasks=tasks,
        ids_path=Path(args.ids),
        out_dir=Path(args.out_dir),
        batch_size=args.batch_size,
        output_format=args.output_format,
        comments=comments,
        cross_references=cross_references,
        ec=ec,
        features=features,
    )
    s.validate()
    return s