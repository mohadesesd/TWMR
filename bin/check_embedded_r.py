#!/usr/bin/env python3
"""Dev helper (v2.2): extract the R script embedded in a Nextflow process (script: \"\"\" ... \"\"\"),
render the Groovy template placeholders with dummy values, and write it out so
that R's parse() can syntax-check it. Usage: check_embedded_r.py file.nf out.R"""
import re, sys

src = open(sys.argv[1]).read()
m = re.search(r'script:\s*\n\s*"""(.*?)"""', src, re.S)
if not m:
    sys.exit("no script block found")
body = m.group(1)

# Groovy GString rendering: \$ -> $, \\ -> \  (order matters: do \\ first)
body = body.replace('\\\\', '\x00BS\x00')
body = body.replace('\\$', '$')
body = body.replace('\x00BS\x00', '\\')

dummies = {
    'params.coloc_rule': 'H4 > 0.8',
    'params.coloc_p12_grid': '1e-5,5e-6,1e-6',
    'params.coloc_min_snps': '50',
    'params.coloc_window_center': 'lead_instrument',
    'params.outcome_type': 'quant',
    "params.outcome_case_prop ?: ''": '',
    'params.locus_plot_set': 'both',
    'params.tissue': 'Whole_Blood',
    'projectDir': '/proj',
}
def render(mo):
    expr = mo.group(1).strip()
    if expr in dummies:
        return dummies[expr]
    # numeric-looking params and val inputs -> 0.8 ; path/val names -> a dummy string
    if any(expr.endswith(k) for k in ('threshold', 'p12', 'window_size', 'cond_threshold', 'power_min', 'h3_strong', 'fdr_target', 'ld_check_r2', 'sample_size', 'fstat_threshold')):
        return '0.8'
    return 'DUMMY_' + re.sub(r'[^A-Za-z0-9_]', '_', expr)
body = re.sub(r'\$\{([^}]*)\}', render, body)
open(sys.argv[2], 'w').write(body)
print(f"wrote {sys.argv[2]} ({body.count(chr(10))} lines)")
