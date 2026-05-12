#!/usr/bin/env python3

from datetime import datetime
from pathlib import Path
import base64
import html
import mimetypes
import re


MARKER_DIR = Path("results/readme_figures/markers")
OUTPUT_HTML = MARKER_DIR / "marker_review.html"
IMAGE_EXTS = {".png", ".jpg", ".jpeg", ".webp", ".gif"}

MARKER_GROUPS = {
    "Rod": ["Rho", "Gnat1", "Pde6a", "Pde6b", "Cnga1", "Nr2e3", "Nrl", "Prph2", "Rom1"],
    "Cone": ["Pde6h", "Arr3", "Opn1mw", "Gnat2", "Pde6c", "Opn1sw"],
    "BC": ["Vsx1", "Vsx2", "Car10", "Prkca", "Sebox", "Scgn", "Cabp5", "Grm6", "Otx2os1", "Gng13", "Nrxn3", "Gabrb3", "Trnp1", "Kcnma1", "Frmd3", "Gm4792", "Nyap2"],
    "AC": ["Tfap2a", "Tfap2b", "Pax6", "Frmd5", "Nrg3", "Elavl3", "Slc32a1", "Slc6a9", "Dcx", "Chat", "Gad1", "Gad2"],
    "MG": ["Rlbp1", "Glul", "Dkk3", "Adamtsl1", "Gpr37", "Abca8a", "Rgs6", "Spc25", "Clu", "Apoe", "Aqp4"],
    "HC": ["Calb1", "Onecut1", "Slc4a3", "Onecut2", "Gm45459", "C1ql1"],
    "RGC": ["Nefl", "Stmn2", "Nrn1", "Pou4f1", "Pou4f2", "Sncg", "Rbpms"],
    "RPE": ["Ttr", "Rdh5", "Rpe65", "Rgr", "Slc16a8"],
    "Astrocyte": ["S100b", "Gfap", "Pax2", "Pdgfra", "Mlc1", "Prdx6"],
    "Microglia": ["Ctss", "C1qa", "Hexb", "Trem2"],
    "Endothelial": ["Cldn5", "Flt1", "Pecam1", "Ly6c1", "Ptprb"],
    "Pericyte": ["Kcnj8", "Acta2", "Pdgfra"],
    "Proliferation": ["Cdk1", "Mki67", "Top2a", "Pcna"],
    "Neurogenic / Progenitor": ["Ascl1", "Neurog2", "Insm1", "Atoh7", "Neurod1", "Sox9", "Lhx2", "Sox2", "Otx2", "Crx", "Olig2", "Foxn4"],
    "Stress / QC": ["Fos", "Jun", "Stat3", "Lcn2", "Malat1"],
}


def clean_gene_name(filename):
    stem = Path(filename).stem
    for pattern in (
        r"^UMAP_Combined_",
        r"^UMAP_RNA_",
        r"^UMAP_ATAC_",
        r"^GeneExpression_",
        r"^GeneScore_",
        r"^marker_",
        r"_expression$",
        r"_GeneExpression$",
        r"_GeneScore$",
        r"_umap$",
        r"_UMAP$",
    ):
        stem = re.sub(pattern, "", stem)
    return stem


def make_anchor(text):
    return re.sub(r"[^A-Za-z0-9_-]+", "_", text)


def build_gene_index(gene_to_file):
    used_genes = set()
    blocks = []

    for group_name, genes in MARKER_GROUPS.items():
        links = []
        for gene in genes:
            used_genes.add(gene)
            label = html.escape(gene)
            if gene in gene_to_file:
                links.append(f'<a href="#{make_anchor(gene)}" data-gene="{label}">{label}</a>')
            else:
                links.append(f'<span class="missing">{label}</span>')

        blocks.append(
            f"""
            <section class="marker-group">
              <h2>{html.escape(group_name)}</h2>
              <div class="gene-links">{" ".join(links)}</div>
            </section>
            """
        )

    extra_genes = sorted(set(gene_to_file) - used_genes, key=str.lower)
    if extra_genes:
        links = [
            f'<a href="#{make_anchor(gene)}" data-gene="{html.escape(gene)}">{html.escape(gene)}</a>'
            for gene in extra_genes
        ]
        blocks.append(
            f"""
            <section class="marker-group">
              <h2>Other</h2>
              <div class="gene-links">{" ".join(links)}</div>
            </section>
            """
        )

    return "\n".join(blocks)


def ordered_genes(gene_to_file):
    ordered = []
    for genes in MARKER_GROUPS.values():
        for gene in genes:
            if gene in gene_to_file and gene not in ordered:
                ordered.append(gene)
    ordered.extend(sorted(set(gene_to_file) - set(ordered), key=str.lower))
    return ordered


def build_cards(gene_to_file):
    cards = []
    for gene in ordered_genes(gene_to_file):
        image_path = gene_to_file[gene]
        safe_gene = html.escape(gene)
        safe_name = html.escape(image_path.name)
        anchor = make_anchor(gene)
        mime_type = mimetypes.guess_type(image_path.name)[0] or "image/png"
        encoded_image = base64.b64encode(image_path.read_bytes()).decode("ascii")
        src = f"data:{mime_type};base64,{encoded_image}"
        cards.append(
            f"""
            <article class="plot-card" id="{anchor}">
              <h2>{safe_gene}</h2>
              <img src="{src}" alt="{safe_gene}" loading="lazy">
              <div class="plot-meta">
                <span>{safe_name}</span>
                <span>embedded</span>
              </div>
            </article>
            """
        )
    return "\n".join(cards)


def main():
    if not MARKER_DIR.exists():
        raise FileNotFoundError(f"Marker directory does not exist: {MARKER_DIR}")

    files = sorted(
        p for p in MARKER_DIR.iterdir()
        if p.is_file() and p.suffix.lower() in IMAGE_EXTS
    )
    if not files:
        raise RuntimeError(f"No marker images found in: {MARKER_DIR}")

    gene_to_file = {clean_gene_name(p.name): p for p in files}
    gene_index = build_gene_index(gene_to_file)
    cards = build_cards(gene_to_file)

    content = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>RNA Marker UMAP Review</title>
  <style>
    :root {{
      --bg: #f4f5f7;
      --panel: #ffffff;
      --text: #1f2933;
      --muted: #667085;
      --line: #d9dee7;
      --link: #155db2;
      --highlight: #ffd34d;
      --highlight-bg: #fff7d1;
    }}

    html {{
      scroll-behavior: smooth;
    }}

    body {{
      margin: 0;
      font-family: Arial, Helvetica, sans-serif;
      background: var(--bg);
      color: var(--text);
    }}

    header {{
      position: sticky;
      top: 0;
      z-index: 10;
      background: rgba(244, 245, 247, 0.96);
      border-bottom: 1px solid var(--line);
      padding: 12px 18px;
      backdrop-filter: blur(6px);
    }}

    h1 {{
      margin: 0 0 4px 0;
      font-size: 22px;
    }}

    .note {{
      margin: 0;
      color: var(--muted);
      font-size: 13px;
    }}

    .layout {{
      display: grid;
      grid-template-columns: minmax(240px, 320px) minmax(0, 1fr);
      gap: 14px;
      padding: 14px;
    }}

    .index-panel {{
      position: sticky;
      top: 78px;
      align-self: start;
      max-height: calc(100vh - 96px);
      overflow: auto;
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 10px;
    }}

    .index-panel h2 {{
      margin: 0 0 8px 0;
      font-size: 15px;
    }}

    .marker-group {{
      border-top: 1px solid #edf0f5;
      padding: 8px 0;
    }}

    .marker-group:first-child {{
      border-top: 0;
      padding-top: 0;
    }}

    .marker-group h2 {{
      margin: 0 0 5px 0;
      font-size: 12px;
      color: #344054;
      text-transform: uppercase;
      letter-spacing: 0.03em;
    }}

    .gene-links {{
      line-height: 1.35;
    }}

    .gene-links a,
    .missing {{
      display: inline-block;
      margin: 0 5px 4px 0;
      font-size: 12px;
    }}

    .gene-links a {{
      color: var(--link);
      text-decoration: none;
      border-bottom: 1px solid transparent;
    }}

    .gene-links a:hover {{
      border-bottom-color: var(--link);
    }}

    .missing {{
      color: #a7acb5;
      text-decoration: line-through;
    }}

    .plot-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fill, minmax(260px, 1fr));
      gap: 10px;
      align-items: start;
    }}

    .plot-card {{
      scroll-margin-top: 90px;
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 8px;
      min-width: 0;
    }}

    .plot-card h2 {{
      display: inline-block;
      margin: 0 0 6px 0;
      padding: 2px 5px;
      border-radius: 5px;
      font-size: 15px;
      line-height: 1.15;
    }}

    .plot-card:target h2 {{
      background: var(--highlight);
      color: #111827;
      box-shadow: 0 0 0 4px var(--highlight-bg);
    }}

    .plot-card img {{
      display: block;
      width: 100%;
      height: auto;
      border: 1px solid #edf0f5;
      border-radius: 6px;
      background: white;
    }}

    .plot-meta {{
      display: flex;
      justify-content: space-between;
      gap: 8px;
      margin-top: 6px;
      color: var(--muted);
      font-size: 11px;
      word-break: break-all;
    }}

    .plot-meta a {{
      color: var(--link);
      text-decoration: none;
      white-space: nowrap;
    }}

    .top-link {{
      position: fixed;
      right: 16px;
      bottom: 16px;
      z-index: 20;
      background: #202938;
      color: white;
      padding: 8px 11px;
      border-radius: 999px;
      text-decoration: none;
      font-size: 12px;
    }}

    @media (max-width: 900px) {{
      .layout {{
        display: block;
      }}

      .index-panel {{
        position: static;
        max-height: none;
        margin-bottom: 12px;
      }}
    }}
  </style>
</head>
<body id="top">
  <header>
    <h1>RNA Marker UMAP Review</h1>
    <p class="note">Generated {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}. Click a gene name to jump to its plot; the selected gene title is highlighted. Images are embedded, so this HTML can be shared as a single file.</p>
  </header>

  <a class="top-link" href="#top">Top</a>

  <main class="layout">
    <aside class="index-panel">
      <h2>Gene Index</h2>
      {gene_index}
    </aside>

    <section class="plot-grid">
      {cards}
    </section>
  </main>
</body>
</html>
"""

    OUTPUT_HTML.write_text(content, encoding="utf-8")
    print(f"Wrote marker review HTML: {OUTPUT_HTML}")


if __name__ == "__main__":
    main()
