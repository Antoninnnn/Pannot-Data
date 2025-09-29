import re
import json
from typing import Dict
import uniprot_utils

_PLACEHOLDER_RE = re.compile(r"\{([A-Za-z0-9_]+)\}")

def _extract_placeholders(template: Dict[str, str]) -> list[str]:
    """
    Extract unique placeholder names from a template.

    Parameters
    ----------
    template : Dict[str, str]
        Template dictionary expected to contain string values under the keys
        'instruction', 'input', and 'output'. Placeholders inside these strings
        have the form `{name}` where `name` matches `[A-Za-z0-9_]+`.

    Returns
    -------
    list[str]
        Ordered list of unique placeholder names (without braces) in the order
        they are first encountered across 'instruction', 'input', and 'output'.

    Notes
    -----
    - Only the three top-level keys ('instruction', 'input', 'output') are scanned.
    - Duplicates are removed while preserving first-seen order.
    """
    seen = []
    for key in ("instruction", "input", "output"):
        if key in template and isinstance(template[key], str):
            for m in _PLACEHOLDER_RE.finditer(template[key]):
                name = m.group(1)
                if name not in seen:
                    seen.append(name)
    return seen

def fill_template_with_uniprot_fields(template: Dict[str, str], protein_id: str) -> Dict[str, str]:
    """
    Fill a template by fetching UniProtKB fields and replacing placeholders.

    Parameters
    ----------
    template : Dict[str, str]
        A dictionary with keys 'instruction', 'input', and 'output'. Each value
        is a string that may contain one or more placeholders of the form `{field}`.
    protein_id : str
        UniProtKB accession ID (e.g., "P69905").

    Returns
    -------
    Dict[str, str]
        A new dictionary with the same keys where every `{placeholder}` in any
        of the three strings is replaced by the **same** raw JSON payload string
        returned from `uniprot_utils.get_protein_fields(protein_id, fields)`.

    Behavior
    --------
    1) Collect the union of all placeholders present in the template.
    2) Join them into a single comma-separated `fields` argument.
    3) Fetch once from UniProtKB and stringify the full JSON response.
    4) Substitute that single string for every `{placeholder}` occurrence.

    Important
    ---------
    - This function deliberately performs **no** mapping between UniProt return
      field names and JSON keys. Every placeholder is replaced with the identical
      raw JSON string for the request.
    - If no placeholders are found, no request is made and the template is returned
      unchanged (aside from creating a new dict instance).
    """
    # 1) Find placeholders and build fields argument
    needed_fields = _extract_placeholders(template)
    fields_arg = ",".join(needed_fields) if needed_fields else ""

    # 2) Fetch once; stringify using default compact JSON format
    data = uniprot_utils.get_protein_fields(protein_id, fields_arg) if fields_arg else {}
    data_str = json.dumps(data, ensure_ascii=False)

    # 3) Replace every {placeholder} with the raw payload string
    def _fill(text: str) -> str:
        if not isinstance(text, str):
            return text
        return _PLACEHOLDER_RE.sub(lambda m: data_str, text)

    return {
        "instruction": _fill(template.get("instruction", "")),
        "input": _fill(template.get("input", "")),
        "output": _fill(template.get("output", "")),
    }