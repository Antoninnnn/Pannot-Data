import time

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

    want_comments = _want(tasks_set, "comments")
    want_xrefs = _want(tasks_set, "cross_references")
    want_ec = _want(tasks_set, "ec")
    want_features = _want(tasks_set, "features")

    # --- load templates ---
    templates = {}

    if want_comments:
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

    if want_xrefs:
        t = s.cross_references
        if t is None:
            raise ValueError("tasks includes 'cross_references' but Settings.cross_references is None")

        templates["go_bp"] = load_templates(t.go_bp)
        templates["go_mf"] = load_templates(t.go_mf)
        templates["go_cc"] = load_templates(t.go_cc)

    if want_ec:
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

    if want_features:
        t = s.features
        if t is None:
            raise ValueError("tasks includes 'features' but Settings.features is None")

        templates["rg_a"] = load_templates(t.rg_a)
        templates["rg_b"] = load_templates(t.rg_b)
        templates["st_a"] = load_templates(t.st_a)
        templates["st_b"] = load_templates(t.st_b)
        templates["st_c"] = load_templates(t.st_c)
        templates["st_d"] = load_templates(t.st_d)

    # --- initialize buffers ---
    BATCH_SIZE = s.batch_size
    protein_ids_batch: list[str] = []

    filled_proteins = 0

    buf_ca, buf_cof, buf_sl, buf_pw = [], [], [], []
    total_samples_ca = total_samples_cof = total_samples_sl = total_samples_pw = 0

    buf_bp, buf_mf, buf_cc = [], [], []
    total_bp = total_mf = total_cc = 0

    buf_ec = []
    total_ec = 0

    buf_rg, buf_st = [], []
    total_rg = total_st = 0

    # --- open files ---
    f_ids = open(s.ids_path, "r", encoding="utf-8")     # input

    out_ctx = []    # outputs (selected tasks only)

    if want_comments:
        out_ctx.extend([
            open(s.ca_out, "a", encoding="utf-8"),
            open(s.cf_out, "a", encoding="utf-8"),
            open(s.sl_out, "a", encoding="utf-8"),
            open(s.pw_out, "a", encoding="utf-8"),
        ])

    if want_xrefs:
        out_ctx.extend([
            open(s.go_bp_out, "a", encoding="utf-8"),
            open(s.go_mf_out, "a", encoding="utf-8"),
            open(s.go_cc_out, "a", encoding="utf-8"),
        ])

    if want_ec:
        out_ctx.append(open(s.ec_out, "a", encoding="utf-8"))

    if want_features:
        out_ctx.extend([
            open(s.rg_out, "a", encoding="utf-8"),
            open(s.st_out, "a", encoding="utf-8"),
        ])

    # --- map file handles ---
    idx = 0
    f_ca = f_cof = f_sl = f_pw = None
    f_bp = f_mf = f_cc = None
    f_ec = None
    f_rg = f_st = None

    if want_comments:
        f_ca, f_cof, f_sl, f_pw = out_ctx[idx:idx+4]
        idx += 4
    if want_xrefs:
        f_bp, f_mf, f_cc = out_ctx[idx:idx+3]
        idx += 3
    if want_ec:
        f_ec = out_ctx[idx]
        idx += 1
    if want_features:
        f_rg, f_st = out_ctx[idx:idx+2]
        idx += 2

    # --- generate samples ---
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

                if want_comments:
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

                if want_xrefs:
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

                if want_ec:
                    ec_samples = gen_ec_samples(
                        meta,
                        templates["ec_a"], templates["ec_b"], templates["ec_c"],
                        templates["ec_d"], templates["ec_e"], templates["ec_f"], templates["ec_g"],
                        templates["ec_h"], templates["ec_i"], templates["ec_j"], templates["ec_k"],
                    )
                    if ec_samples:
                        buf_ec.extend(ec_samples); total_ec += len(ec_samples)

                if want_features:
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

            # --- flush batch ---
            if want_comments:
                write_samples(f_ca, buf_ca, s.output_format); buf_ca.clear()
                write_samples(f_cof, buf_cof, s.output_format); buf_cof.clear()
                write_samples(f_sl, buf_sl, s.output_format); buf_sl.clear()
                write_samples(f_pw, buf_pw, s.output_format); buf_pw.clear()

            if want_xrefs:
                write_samples(f_bp, buf_bp, s.output_format); buf_bp.clear()
                write_samples(f_mf, buf_mf, s.output_format); buf_mf.clear()
                write_samples(f_cc, buf_cc, s.output_format); buf_cc.clear()

            if want_ec:
                write_samples(f_ec, buf_ec, s.output_format); buf_ec.clear()

            if want_features:
                write_samples(f_rg, buf_rg, s.output_format); buf_rg.clear()
                write_samples(f_st, buf_st, s.output_format); buf_st.clear()

            # --- print progress ---
            parts = [f"Filled with {filled_proteins} proteins"]
            if want_comments:
                total_c = total_samples_ca + total_samples_cof + total_samples_sl + total_samples_pw
                parts.append(f"comments={total_c} (CA={total_samples_ca}, CF={total_samples_cof}, SL={total_samples_sl}, PW={total_samples_pw})")
            if want_xrefs:
                total_go = total_bp + total_mf + total_cc
                parts.append(f"go={total_go} (BP={total_bp}, MF={total_mf}, CC={total_cc})")
            if want_ec:
                parts.append(f"ec={total_ec}")
            if want_features:
                total_f = total_rg + total_st
                parts.append(f"features={total_f} (RG={total_rg}, ST={total_st})")

            print(" | ".join(parts))

            protein_ids_batch.clear()

        # --- handle last batch ---
        if protein_ids_batch:
            entries = get_proteins_info(protein_ids_batch)

            for data in entries:
                meta = extract_uniprot_meta(data)

                if want_comments:
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

                if want_xrefs:
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

                if want_ec:
                    ec_samples = gen_ec_samples(
                        meta,
                        templates["ec_a"], templates["ec_b"], templates["ec_c"], templates["ec_d"],
                        templates["ec_e"], templates["ec_f"], templates["ec_g"], templates["ec_h"],
                        templates["ec_i"], templates["ec_j"], templates["ec_k"],
                    )
                    if ec_samples:
                        buf_ec.extend(ec_samples); total_ec += len(ec_samples)

                if want_features:
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

            # --- flush last batch ---
            if want_comments:
                write_samples(f_ca, buf_ca, s.output_format); buf_ca.clear()
                write_samples(f_cof, buf_cof, s.output_format); buf_cof.clear()
                write_samples(f_sl, buf_sl, s.output_format); buf_sl.clear()
                write_samples(f_pw, buf_pw, s.output_format); buf_pw.clear()

            if want_xrefs:
                write_samples(f_bp, buf_bp, s.output_format); buf_bp.clear()
                write_samples(f_mf, buf_mf, s.output_format); buf_mf.clear()
                write_samples(f_cc, buf_cc, s.output_format); buf_cc.clear()

            if want_ec:
                write_samples(f_ec, buf_ec, s.output_format); buf_ec.clear()

            if want_features:
                write_samples(f_rg, buf_rg, s.output_format); buf_rg.clear()
                write_samples(f_st, buf_st, s.output_format); buf_st.clear()

    finally:
        for f in out_ctx:
            try:
                f.close()
            except Exception:
                pass

    # --- summary ---
    parts = [f"Done. Filled with {filled_proteins} proteins"]
    if want_comments:
        total_c = total_samples_ca + total_samples_cof + total_samples_sl + total_samples_pw
        parts.append(f"comments={total_c} (CA={total_samples_ca}, CF={total_samples_cof}, SL={total_samples_sl}, PW={total_samples_pw})")
    if want_xrefs:
        total_go = total_bp + total_mf + total_cc
        parts.append(f"go={total_go} (BP={total_bp}, MF={total_mf}, CC={total_cc})")
    if want_ec:
        parts.append(f"ec={total_ec}")
    if want_features:
        total_f = total_rg + total_st
        parts.append(f"features={total_f} (RG={total_rg}, ST={total_st})")

    print(" | ".join(parts))
    print(f"TOTAL TIME: {time.time() - t_start:.2f}s")