import sys
import time
from pathlib import Path

# Add src/ to Python path (project root–relative, safe)
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_PATH = PROJECT_ROOT / "src"
sys.path.insert(0, str(SRC_PATH))

from Composer.config.cli import parse_settings

from Composer.uniprot_utils import (
    get_proteins_info,
    extract_uniprot_meta,
    extract_uniprot_comments,
    extract_uniprot_cross_references,
    extract_uniprot_features,
)
from Composer.template_utils import (
    load_templates,
    write_samples,
    gen_catalytic_activity_samples,
    gen_cofactor_samples,
    gen_subcellular_location_samples,
    gen_pathway_samples,
    gen_go_samples,
    gen_ec_samples,
    gen_regions_samples,
    gen_sites_samples,
)


def _want(tasks_set: set[str], name: str) -> bool:
    return name in tasks_set


if __name__ == "__main__":
    t_start = time.time()
    s = parse_settings()

    tasks = tuple(s.tasks)
    tasks_set = set(tasks)

    # ----------------------
    # Load templates (only for selected tasks)
    # ----------------------
    templates = {}

    if _want(tasks_set, "comments"):
        t = s.comments
        if t is None:
            raise ValueError("tasks includes 'comments' but Settings.comments is None")

        templates["ca"] = load_templates(t.ca_a)
        templates["cof_a"] = load_templates(t.cf_a)
        templates["cof_b"] = load_templates(t.cf_b)
        templates["sl_a"] = load_templates(t.sl_a)
        templates["sl_b"] = load_templates(t.sl_b)
        templates["sl_c"] = load_templates(t.sl_c)
        templates["sl_d"] = load_templates(t.sl_d)
        templates["pw_a"] = load_templates(t.pw_a)
        templates["pw_b"] = load_templates(t.pw_b)
        templates["pw_c"] = load_templates(t.pw_c)
        templates["pw_d"] = load_templates(t.pw_d)

    if _want(tasks_set, "cross_references"):
        t = s.cross_references
        if t is None:
            raise ValueError("tasks includes 'cross_references' but Settings.cross_references is None")

        templates["go_bp"] = load_templates(t.go_bp)
        templates["go_mf"] = load_templates(t.go_mf)
        templates["go_cc"] = load_templates(t.go_cc)

    if _want(tasks_set, "ec"):
        t = s.ec
        if t is None:
            raise ValueError("tasks includes 'ec' but Settings.ec is None")

        templates["ec_a"] = load_templates(t.ec_a)
        templates["ec_b"] = load_templates(t.ec_b)
        templates["ec_c"] = load_templates(t.ec_c)
        templates["ec_d"] = load_templates(t.ec_d)
        templates["ec_e"] = load_templates(t.ec_e)
        templates["ec_f"] = load_templates(t.ec_f)
        templates["ec_g"] = load_templates(t.ec_g)
        templates["ec_h"] = load_templates(t.ec_h)
        templates["ec_i"] = load_templates(t.ec_i)
        templates["ec_j"] = load_templates(t.ec_j)
        templates["ec_k"] = load_templates(t.ec_k)

    if _want(tasks_set, "features"):
        t = s.features
        if t is None:
            raise ValueError("tasks includes 'features' but Settings.features is None")

        templates["rg_a"] = load_templates(t.rg_a)
        templates["rg_b"] = load_templates(t.rg_b)
        templates["st_a"] = load_templates(t.st_a)
        templates["st_b"] = load_templates(t.st_b)
        templates["st_c"] = load_templates(t.st_c)
        templates["st_d"] = load_templates(t.st_d)

    # ----------------------
    # Buffers + counters
    # ----------------------
    BATCH_SIZE = s.batch_size
    protein_ids_batch: list[str] = []

    filled_proteins = 0

    # comments
    buf_ca, buf_cof, buf_sl, buf_pw = [], [], [], []
    total_samples_ca = total_samples_cof = total_samples_sl = total_samples_pw = 0

    # cross_references
    buf_bp, buf_mf, buf_cc = [], [], []
    total_bp = total_mf = total_cc = 0

    # ec
    buf_ec = []
    total_ec = 0

    # features
    buf_rg, buf_st = [], []
    total_rg = total_st = 0

    # ----------------------
    # Open output files only for selected tasks
    # ----------------------
    out_ctx = []
    if _want(tasks_set, "comments"):
        out_ctx.extend([
            open(s.ca_out, "a", encoding="utf-8"),
            open(s.cf_out, "a", encoding="utf-8"),
            open(s.sl_out, "a", encoding="utf-8"),
            open(s.pw_out, "a", encoding="utf-8"),
        ])
    if _want(tasks_set, "cross_references"):
        out_ctx.extend([
            open(s.go_bp_out, "a", encoding="utf-8"),
            open(s.go_mf_out, "a", encoding="utf-8"),
            open(s.go_cc_out, "a", encoding="utf-8"),
        ])
    if _want(tasks_set, "ec"):
        out_ctx.append(open(s.ec_out, "a", encoding="utf-8"))
    if _want(tasks_set, "features"):
        out_ctx.extend([
            open(s.rg_out, "a", encoding="utf-8"),
            open(s.st_out, "a", encoding="utf-8"),
        ])

    # ids file always needed
    out_ctx.append(open(s.ids_path, "r", encoding="utf-8"))

    # Map handles by task (kept simple, positional based on how we append above)
    idx = 0
    f_ca = f_cof = f_sl = f_pw = None
    f_bp = f_mf = f_cc = None
    f_ec = None
    f_rg = f_st = None

    if _want(tasks_set, "comments"):
        f_ca, f_cof, f_sl, f_pw = out_ctx[idx:idx+4]
        idx += 4
    if _want(tasks_set, "cross_references"):
        f_bp, f_mf, f_cc = out_ctx[idx:idx+3]
        idx += 3
    if _want(tasks_set, "ec"):
        f_ec = out_ctx[idx]
        idx += 1
    if _want(tasks_set, "features"):
        f_rg, f_st = out_ctx[idx:idx+2]
        idx += 2

    f_ids = out_ctx[-1]

    try:
        for line in f_ids:
            protein_id = line.strip()
            if not protein_id or protein_id.startswith("#"):
                continue

            protein_ids_batch.append(protein_id)

            if len(protein_ids_batch) < BATCH_SIZE:
                continue

            entries = get_proteins_info(protein_ids_batch)

            for data in entries:
                meta = extract_uniprot_meta(data)

                if _want(tasks_set, "comments"):
                    comments = extract_uniprot_comments(data)

                    ca_samples = gen_catalytic_activity_samples(comments, meta, templates["ca"])
                    cof_samples = gen_cofactor_samples(comments, meta, templates["cof_a"], templates["cof_b"])
                    sl_samples = gen_subcellular_location_samples(
                        comments, meta,
                        templates["sl_a"], templates["sl_b"], templates["sl_c"], templates["sl_d"]
                    )
                    pw_samples = gen_pathway_samples(
                        comments, meta,
                        templates["pw_a"], templates["pw_b"], templates["pw_c"], templates["pw_d"]
                    )

                    if ca_samples:
                        buf_ca.extend(ca_samples); total_samples_ca += len(ca_samples)
                    if cof_samples:
                        buf_cof.extend(cof_samples); total_samples_cof += len(cof_samples)
                    if sl_samples:
                        buf_sl.extend(sl_samples); total_samples_sl += len(sl_samples)
                    if pw_samples:
                        buf_pw.extend(pw_samples); total_samples_pw += len(pw_samples)

                if _want(tasks_set, "cross_references"):
                    cross_refs = extract_uniprot_cross_references(data)

                    go_samples = gen_go_samples(
                        cross_refs,
                        meta,
                        templates["go_bp"], templates["go_mf"], templates["go_cc"],
                    )

                    if go_samples.get("BP"):
                        buf_bp.append(go_samples["BP"]); total_bp += 1
                    if go_samples.get("MF"):
                        buf_mf.append(go_samples["MF"]); total_mf += 1
                    if go_samples.get("CC"):
                        buf_cc.append(go_samples["CC"]); total_cc += 1

                if _want(tasks_set, "ec"):
                    ec_samples = gen_ec_samples(
                        meta,
                        templates["ec_a"], templates["ec_b"], templates["ec_c"],
                        templates["ec_d"], templates["ec_e"], templates["ec_f"], templates["ec_g"],
                        templates["ec_h"], templates["ec_i"], templates["ec_j"], templates["ec_k"],
                    )
                    if ec_samples:
                        buf_ec.extend(ec_samples); total_ec += len(ec_samples)

                if _want(tasks_set, "features"):
                    features = extract_uniprot_features(data)

                    rg_samples = gen_regions_samples(features, meta, templates["rg_a"], templates["rg_b"])
                    st_samples = gen_sites_samples(
                        features, meta,
                        templates["st_a"], templates["st_b"], templates["st_c"], templates["st_d"]
                    )

                    if rg_samples:
                        buf_rg.extend(rg_samples); total_rg += len(rg_samples)
                    if st_samples:
                        buf_st.extend(st_samples); total_st += len(st_samples)

                filled_proteins += 1

            # Flush once per batch (only for enabled tasks)
            if _want(tasks_set, "comments"):
                write_samples(f_ca, buf_ca, s.output_format); buf_ca.clear()
                write_samples(f_cof, buf_cof, s.output_format); buf_cof.clear()
                write_samples(f_sl, buf_sl, s.output_format); buf_sl.clear()
                write_samples(f_pw, buf_pw, s.output_format); buf_pw.clear()

            if _want(tasks_set, "cross_references"):
                write_samples(f_bp, buf_bp, s.output_format); buf_bp.clear()
                write_samples(f_mf, buf_mf, s.output_format); buf_mf.clear()
                write_samples(f_cc, buf_cc, s.output_format); buf_cc.clear()

            if _want(tasks_set, "ec"):
                write_samples(f_ec, buf_ec, s.output_format); buf_ec.clear()

            if _want(tasks_set, "features"):
                write_samples(f_rg, buf_rg, s.output_format); buf_rg.clear()
                write_samples(f_st, buf_st, s.output_format); buf_st.clear()

            # Print progress
            parts = [f"Filled with {filled_proteins} proteins"]
            if _want(tasks_set, "comments"):
                total_c = total_samples_ca + total_samples_cof + total_samples_sl + total_samples_pw
                parts.append(f"comments={total_c} (CA={total_samples_ca}, CF={total_samples_cof}, SL={total_samples_sl}, PW={total_samples_pw})")
            if _want(tasks_set, "cross_references"):
                total_go = total_bp + total_mf + total_cc
                parts.append(f"go={total_go} (BP={total_bp}, MF={total_mf}, CC={total_cc})")
            if _want(tasks_set, "ec"):
                parts.append(f"ec={total_ec}")
            if _want(tasks_set, "features"):
                total_f = total_rg + total_st
                parts.append(f"features={total_f} (RG={total_rg}, ST={total_st})")

            print(" | ".join(parts))

            protein_ids_batch.clear()

        # Handle last batch (<BATCH_SIZE)
        if protein_ids_batch:
            entries = get_proteins_info(protein_ids_batch)

            for data in entries:
                meta = extract_uniprot_meta(data)

                if _want(tasks_set, "comments"):
                    comments = extract_uniprot_comments(data)

                    ca_samples = gen_catalytic_activity_samples(comments, meta, templates["ca"])
                    cof_samples = gen_cofactor_samples(comments, meta, templates["cof_a"], templates["cof_b"])
                    sl_samples = gen_subcellular_location_samples(
                        comments, meta,
                        templates["sl_a"], templates["sl_b"], templates["sl_c"], templates["sl_d"]
                    )
                    pw_samples = gen_pathway_samples(
                        comments, meta,
                        templates["pw_a"], templates["pw_b"], templates["pw_c"], templates["pw_d"]
                    )

                    if ca_samples:
                        buf_ca.extend(ca_samples); total_samples_ca += len(ca_samples)
                    if cof_samples:
                        buf_cof.extend(cof_samples); total_samples_cof += len(cof_samples)
                    if sl_samples:
                        buf_sl.extend(sl_samples); total_samples_sl += len(sl_samples)
                    if pw_samples:
                        buf_pw.extend(pw_samples); total_samples_pw += len(pw_samples)

                if _want(tasks_set, "cross_references"):
                    cross_refs = extract_uniprot_cross_references(data)

                    go_samples = gen_go_samples(
                        cross_refs,
                        meta,
                        templates["go_bp"], templates["go_mf"], templates["go_cc"],
                    )

                    if go_samples.get("BP"):
                        buf_bp.append(go_samples["BP"]); total_bp += 1
                    if go_samples.get("MF"):
                        buf_mf.append(go_samples["MF"]); total_mf += 1
                    if go_samples.get("CC"):
                        buf_cc.append(go_samples["CC"]); total_cc += 1

                if _want(tasks_set, "ec"):
                    ec_samples = gen_ec_samples(
                        meta,
                        templates["ec_a"], templates["ec_b"], templates["ec_c"],
                        templates["ec_d"], templates["ec_e"], templates["ec_f"], templates["ec_g"],
                        templates["ec_h"], templates["ec_i"], templates["ec_j"], templates["ec_k"],
                    )
                    if ec_samples:
                        buf_ec.extend(ec_samples); total_ec += len(ec_samples)

                if _want(tasks_set, "features"):
                    features = extract_uniprot_features(data)

                    rg_samples = gen_regions_samples(features, meta, templates["rg_a"], templates["rg_b"])
                    st_samples = gen_sites_samples(
                        features, meta,
                        templates["st_a"], templates["st_b"], templates["st_c"], templates["st_d"]
                    )

                    if rg_samples:
                        buf_rg.extend(rg_samples); total_rg += len(rg_samples)
                    if st_samples:
                        buf_st.extend(st_samples); total_st += len(st_samples)

                filled_proteins += 1

            # Flush last batch
            if _want(tasks_set, "comments"):
                write_samples(f_ca, buf_ca, s.output_format); buf_ca.clear()
                write_samples(f_cof, buf_cof, s.output_format); buf_cof.clear()
                write_samples(f_sl, buf_sl, s.output_format); buf_sl.clear()
                write_samples(f_pw, buf_pw, s.output_format); buf_pw.clear()

            if _want(tasks_set, "cross_references"):
                write_samples(f_bp, buf_bp, s.output_format); buf_bp.clear()
                write_samples(f_mf, buf_mf, s.output_format); buf_mf.clear()
                write_samples(f_cc, buf_cc, s.output_format); buf_cc.clear()

            if _want(tasks_set, "ec"):
                write_samples(f_ec, buf_ec, s.output_format); buf_ec.clear()

            if _want(tasks_set, "features"):
                write_samples(f_rg, buf_rg, s.output_format); buf_rg.clear()
                write_samples(f_st, buf_st, s.output_format); buf_st.clear()

    finally:
        for f in out_ctx:
            try:
                f.close()
            except Exception:
                pass

    # Summary
    parts = [f"Done. Filled with {filled_proteins} proteins"]
    if _want(tasks_set, "comments"):
        total_c = total_samples_ca + total_samples_cof + total_samples_sl + total_samples_pw
        parts.append(f"comments={total_c} (CA={total_samples_ca}, CF={total_samples_cof}, SL={total_samples_sl}, PW={total_samples_pw})")
    if _want(tasks_set, "cross_references"):
        total_go = total_bp + total_mf + total_cc
        parts.append(f"go={total_go} (BP={total_bp}, MF={total_mf}, CC={total_cc})")
    if _want(tasks_set, "ec"):
        parts.append(f"ec={total_ec}")
    if _want(tasks_set, "features"):
        total_f = total_rg + total_st
        parts.append(f"features={total_f} (RG={total_rg}, ST={total_st})")

    print(" | ".join(parts))
    print(f"TOTAL TIME: {time.time() - t_start:.2f}s")