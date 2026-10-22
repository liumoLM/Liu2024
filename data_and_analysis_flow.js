// Build data_and_analysis_flow.pptx: how data flow from unified_indels_data
// through Liu2026_code_and_data to the paper's figures, tables, numbers, and
// indelsig.net. Same visual vocabulary as topography_pipeline.pptx.
// Run with pptxgenjs installed (NODE_PATH pointing at its node_modules):
//   NODE_PATH=<dir>/node_modules node data_and_analysis_flow.js data_and_analysis_flow.pptx
// The theme step uses apply_theme.js from the Claude pptx skill (path below).
const pptxgen = require("pptxgenjs");
const { applyTheme } = require("/home/steve/.claude/skills/synced/7ecff367-0fa1-499d-b6f2-b54ca49672e7_58378af3-09cf-4e76-af92-9fc218b3e24b/pptx/scripts/apply_theme.js");

const THEME = {
  name: "Liu 2026 pipeline",
  headFontFace: "Cambria",
  bodyFontFace: "Calibri",
  colors: {
    dk1: "1F2A30", lt1: "FFFFFF", dk2: "3B4A54", lt2: "EEF2F4",
    accent1: "2E7D4F", accent2: "C0392B", accent3: "1F5F8B",
    accent4: "C98A1B", accent5: "7A8A94", accent6: "D9E7DE",
    hlink: "1F5F8B", folHlink: "6D4C8B",
  },
};
const HEX = THEME.colors;
const STATUS = "status 2026-10-22";

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE"; // 13.33 x 7.5 in
pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
pres.title = "Liu et al. 2026: data and analysis flow";
const C = pres.SchemeColor;

pres.defineSlideMaster({
  title: "FLOW",
  background: { color: HEX.lt1 },
  objects: [
    { text: { text: `Liu et al. 2026 · data and analysis flow · ${STATUS}`,
      options: { x: 0.5, y: 7.05, w: 8, h: 0.3, fontSize: 9, color: HEX.accent5 } } },
  ],
  slideNumber: { x: 12.4, y: 7.05, w: 0.5, h: 0.3, fontSize: 9, color: HEX.accent5 },
});

// ---- box vocabulary (as in topography_pipeline.pptx) ----------------------
const KIND = {
  script:  { fill: HEX.lt1, line: HEX.dk1, text: HEX.dk1, mono: true },
  data:    { fill: HEX.accent1, line: HEX.accent1, text: HEX.lt1 },
  other:   { fill: HEX.accent4, line: HEX.accent4, text: HEX.lt1 },
  missing: { fill: HEX.accent2, line: HEX.accent2, text: HEX.lt1 },
  output:  { fill: HEX.accent3, line: HEX.accent3, text: HEX.lt1 },
  stale:   { fill: HEX.lt1, line: HEX.accent2, text: HEX.accent2, dash: true },
  hand:    { fill: HEX.lt1, line: HEX.accent3, text: HEX.accent3, dash: true },
  nocode:  { fill: HEX.lt1, line: HEX.accent2, text: HEX.accent2, dash: true, mono: true },
};

function box(slide, kind, x, y, w, h, text, opts = {}) {
  const k = KIND[kind];
  slide.addText(text, {
    shape: pres.shapes.ROUNDED_RECTANGLE, rectRadius: 0.06,
    x, y, w, h, isTextBox: true,
    fill: { color: k.fill },
    line: { color: k.line, width: 1, dashType: k.dash ? "dash" : "solid" },
    color: k.text, fontSize: opts.fontSize || 11,
    fontFace: k.mono ? "Courier New" : "Calibri",
    bold: opts.bold || false, align: opts.align || "center", valign: "middle",
    margin: 3, objectName: opts.name || text.toString().slice(0, 40),
  });
}

function arrow(slide, x1, y1, x2, y2) {
  slide.addShape(pres.shapes.LINE, {
    x: Math.min(x1, x2), y: Math.min(y1, y2),
    w: Math.abs(x2 - x1) || 0.001, h: Math.abs(y2 - y1) || 0.001,
    flipV: y2 < y1,
    line: { color: HEX.accent5, width: 1.25, endArrowType: "triangle" },
  });
}

function note(slide, x, y, w, text, color = HEX.accent3) {
  slide.addText(text, { x, y, w, h: 0.3, fontSize: 10, italic: true, color, isTextBox: true, margin: 0 });
}

function legend(slide, kinds) {
  const labels = {
    script: "Script", data: "Data in data repo", other: "Data, other folder",
    missing: "Data missing", output: "Built output", stale: "Changes with new set",
    hand: "Edited by hand", nocode: "Code not in repos",
  };
  let x = 0.5;
  for (const k of kinds) {
    box(slide, k, x, 6.6, 0.35, 0.22, "", { name: `legend ${k}` });
    slide.addText(labels[k], { x: x + 0.42, y: 6.56, w: 1.6, h: 0.3, fontSize: 10, color: HEX.dk1, isTextBox: true, margin: 0 });
    x += 2.05;
  }
}

function newSlide(title, subtitle, section) {
  const s = pres.addSlide({ masterName: "FLOW", sectionTitle: section });
  s.addText(title, { x: 0.5, y: 0.3, w: 12.3, h: 0.6, fontSize: 28, bold: true, color: HEX.dk1,
    fontFace: "Cambria", align: "center", isTextBox: true, margin: 0, objectName: "Title" });
  s.addText(subtitle, { x: 0.5, y: 0.95, w: 12.3, h: 0.45, fontSize: 14, color: HEX.dk2,
    isTextBox: true, margin: 0, objectName: "Subtitle" });
  return s;
}

// One row: inputs -> script -> output, with an optional note at the right.
function flowRow(slide, y, inputs, script, output, opts = {}) {
  const h = opts.h || 0.5;
  const xIn = 0.5, wIn = 3.9, xS = 4.85, wS = 3.3, xO = 8.6, wO = 2.2;
  box(slide, opts.inKind || "data", xIn, y, wIn, h, inputs, { fontSize: 10 });
  arrow(slide, xIn + wIn, y + h / 2, xS, y + h / 2);
  box(slide, opts.scriptKind || "script", xS, y, wS, h, script, { fontSize: 10 });
  arrow(slide, xS + wS, y + h / 2, xO, y + h / 2);
  box(slide, opts.outKind || "output", xO, y, wO, h, output, { fontSize: 11, bold: true });
  if (opts.note) note(slide, xO + wO + 0.15, y + h / 2 - 0.15, 12.85 - (xO + wO + 0.15), opts.note, opts.noteColor);
}

// ---- Slide 1: the big picture ---------------------------------------------
{
  pres.addSection({ title: "Overview" });
  const s = newSlide(
    "From raw data to the paper and indelsig.net",
    "Three repositories. The data repo holds the finalized data, the code repo turns them into every figure and table, and the paper repo assembles them.",
    "Overview");
  const col = (x, w, head) => s.addText(head, { x, y: 1.5, w, h: 0.35, fontSize: 14, bold: true, color: HEX.dk1, isTextBox: true, margin: 0 });
  col(0.5, 2.6, "Upstream (not rerun here)");
  col(3.55, 2.7, "unified_indels_data");
  col(6.75, 2.9, "Liu2026_code_and_data");
  col(10.1, 2.75, "Outputs");

  box(s, "other", 0.5, 1.95, 2.6, 0.7, "HMF and PCAWG VCFs\nmetadata (Hartwig, ICGC)", { fontSize: 10 });
  box(s, "script", 0.5, 2.9, 2.6, 0.7, "mSigHdp, SigProfiler\nExtractor (extraction)", { fontSize: 10 });
  box(s, "script", 0.5, 3.85, 2.6, 0.7, "SigProfilerAssignment\n+ mSigAct (Mo)", { fontSize: 10 });
  box(s, "other", 0.5, 4.8, 2.6, 0.7, "external data\nDegasperi 2022, Reijns 2022", { fontSize: 10 });

  box(s, "data", 3.55, 1.95, 2.7, 0.7, "spectra 83 / 89 / 476\nsample_info.tsv", { fontSize: 10 });
  box(s, "data", 3.55, 2.9, 2.7, 0.7, "signatures 83 / 89 / 476\nconnection_table.tsv", { fontSize: 10 });
  box(s, "data", 3.55, 3.85, 2.7, 0.7, "assignments 83 / 89\n(new set 2026-10-19)", { fontSize: 10 });
  box(s, "data", 3.55, 4.8, 2.7, 0.7, "external_data/\nraw extraction runs", { fontSize: 10 });
  for (const y of [2.3, 3.25, 4.2, 5.15]) arrow(s, 3.1, y, 3.55, y);

  box(s, "script", 6.75, 1.95, 2.9, 1.15, "build_fig1/ … build_fig6/\nbuild_fig_s1/ … build_fig_s6/\nbuild_sup_table_s1/ … s17/", { fontSize: 10 });
  box(s, "script", 6.75, 3.3, 2.9, 0.7, "ID83_ / ID89_topography_\nanalysis_code/ (Mini)", { fontSize: 10 });
  box(s, "script", 6.75, 4.2, 2.9, 0.55, "vignette/", { fontSize: 10 });
  box(s, "script", 6.75, 4.95, 2.9, 0.55, "ms_numbers.R", { fontSize: 10 });
  for (const y of [2.3, 3.25, 4.2, 5.15]) arrow(s, 6.25, y, 6.75, Math.min(Math.max(y, 2.2), 5.2));

  box(s, "output", 10.1, 1.95, 2.75, 0.5, "Figures 1–6", { bold: true });
  box(s, "output", 10.1, 2.55, 2.75, 0.5, "Figures S1–S6", { bold: true });
  box(s, "output", 10.1, 3.15, 2.75, 0.5, "Tables S1–S17", { bold: true });
  box(s, "output", 10.1, 4.2, 2.75, 0.55, "indelsig.net pages", { bold: true });
  box(s, "output", 10.1, 4.95, 2.75, 0.55, "numbers in ms.qmd", { bold: true });
  arrow(s, 9.65, 2.5, 10.1, 2.2); arrow(s, 9.65, 2.55, 10.1, 2.8); arrow(s, 9.65, 2.6, 10.1, 3.4);
  arrow(s, 9.65, 3.65, 10.1, 3.45);
  arrow(s, 9.65, 4.48, 10.1, 4.48); arrow(s, 9.65, 5.23, 10.1, 5.23);

  s.addText("Code finds the data with find_data.R (or UNIFIED_INDELS_DATA). The paper repo copies each built file under the same name, and check_main_figures.R, check_sup_figures.R and check_sup_tables.R flag copies that are out of date.",
    { x: 0.5, y: 5.75, w: 12.3, h: 0.6, fontSize: 12, color: HEX.dk2, isTextBox: true, margin: 0 });
  legend(s, ["script", "data", "other", "output"]);
  s.addNotes("Overview of the three repositories. Upstream steps (variant calling, extraction, assignment) are not rerun by the build scripts. Their results are frozen in unified_indels_data.");
}

// ---- Slide 2: the data repo ------------------------------------------------
{
  pres.addSection({ title: "Data" });
  const s = newSlide(
    "unified_indels_data: what each file is and who made it",
    "Every build script reads these files. When one changes, the figures and tables downstream of it must be rebuilt (slide 7).",
    "Data");
  const rows = [
    ["HMF and PCAWG indel VCFs, clipped to repeat length ≤ 9", "catalog generation", "liu_et_al_{83,89,476}_spectra.tsv", "nocode"],
    ["Hartwig metadata, ICGC donor table, MSIseq calls", "build_all_meta.R (data repo)", "sample_info.tsv"],
    ["spectra (de novo extraction)", "mSigHdp, SigProfilerExtractor", "liu_et_al_{83,89,476}_signatures.tsv"],
    ["signatures + spectra (best-matching tumors)", "make_connection_table.R", "connection_table.tsv"],
    ["172-type spectra, signatures", "SigProfilerAssignment, then mSigAct for 731 tumors", "liu_et_al_{83,89}_assignment.tsv"],
    ["Degasperi 2022 SBS exposures; Reijns 2022 spectra", "downloaded", "external_data/"],
  ];
  let y = 1.6;
  for (const [inp, scr, out, sk] of rows) {
    flowRow(s, y, inp, scr, out, { inKind: inp.startsWith("spectra") || inp.startsWith("signatures") || inp.startsWith("172") ? "data" : "other", outKind: "data", h: 0.6, scriptKind: sk });
    y += 0.78;
  }
  note(s, 8.6, y - 0.05, 4.2, "Mo's assignments of 2026-10-19 are in assignments_2026-10-19/ and are not yet the canonical files.", HEX.accent2);
  legend(s, ["script", "nocode", "data", "other"]);
}

// ---- Slide 3: main figures --------------------------------------------------
{
  pres.addSection({ title: "Figures and tables" });
  const s = newSlide(
    "Main figures",
    "Each build_figN/ writes figN.pdf. The paper's main_figures/ copies it, and main_figures_manifest.tsv records the source.",
    "Figures and tables");
  const h = 0.6, dy = 0.78;
  flowRow(s, 1.6, "signatures 83 / 89 / 476", "build_fig1/build_fig1.R", "Figure 1", { h });
  flowRow(s, 1.6 + dy, "liu_et_al_83_signatures.tsv", "build_fig2/build_fig2.R", "Figure 2", { h });
  flowRow(s, 1.6 + 2 * dy, "liu_et_al_89_signatures.tsv", "build_fig3/build_fig3.R", "Figure 3", { h });
  flowRow(s, 1.6 + 3 * dy, "assignments 83 / 89, sample_info, spectra, Degasperi SBS", "build_fig4/fig4_combined.R\n(enrichment + SBS-ID clustering)", "Figure 4", { h, outKind: "stale", note: "4A loses 3 dots, gains 1" , noteColor: HEX.accent2});
  flowRow(s, 1.6 + 4 * dy, "topography r.04 matrices (Mini)", "build_fig5/ (see topography deck)", "Figure 5", { h, inKind: "data", outKind: "stale", note: "needs Mo's new partial credit", noteColor: HEX.accent2 });
  flowRow(s, 1.6 + 5 * dy, "signatures, Reijns 2022 spectra", "build_fig6/Plot signatures and motifs.qmd", "Figure 6", { h });
  legend(s, ["script", "data", "output", "stale"]);
}

// ---- Slide 4: supplementary figures ----------------------------------------
{
  const s = newSlide(
    "Supplementary figures",
    "Figures S2 and S5 are edited by hand from a _draft source. sup_figures_manifest.tsv maps each one to its source.",
    "Figures and tables");
  const h = 0.6, dy = 0.78;
  flowRow(s, 1.6, "assignments, signatures, spectra, sample_info", "build_fig_s1/Signature Assignment.qmd", "Figure S1", { h, outKind: "stale", note: "colors and dots change", noteColor: HEX.accent2 });
  flowRow(s, 1.6 + dy, "Tables S2, S4, S6 signatures", "build_fig_s2/make_figure_s2_layout.R", "Figure S2", { h, outKind: "hand" });
  flowRow(s, 1.6 + 2 * dy, "topography Fisher tests (r.05, r.06)", "build_fig_s3/r.08_bubbleplot.R", "Figure S3", { h, outKind: "stale", note: "needs new partial credit", noteColor: HEX.accent2 });
  flowRow(s, 1.6 + 3 * dy, "replication-timing trends (r.09)", "build_fig_s4/r.09_barplot_final.R", "Figure S4", { h, outKind: "stale", note: "needs new partial credit", noteColor: HEX.accent2 });
  flowRow(s, 1.6 + 4 * dy, "Reijns 2022 spectra", "build_fig_s5/make_figure_s5_layout.R", "Figure S5", { h, outKind: "hand" });
  flowRow(s, 1.6 + 5 * dy, "assignments, sample_info", "build_fig_s6/build_fig_s6_stats.R\n+ make_figure_s6_layout.R", "Figure S6", { h, outKind: "stale", note: "sex hits change", noteColor: HEX.accent2 });
  legend(s, ["script", "data", "output", "stale", "hand"]);
}

// ---- Slide 5: supplementary tables -----------------------------------------
{
  const s = newSlide(
    "Supplementary tables",
    "sup_tables_manifest.tsv maps each table to its source, and check_sup_tables.R compares them cell by cell.",
    "Figures and tables");
  const h = 0.5, dy = 0.6;
  let y = 1.5;
  const R = (inp, scr, out, o = {}) => { flowRow(s, y, inp, scr, out, Object.assign({ h }, o)); y += dy; };
  R("HMF and PCAWG annotated VCFs", "build_sup_table_s1/ (3 steps)", "Table S1", { inKind: "other" });
  R("raw extraction runs, signatures", "build_sup_table_s2/ … s8/", "Tables S2–S8");
  R("89 assignment + spectra, coding variants", "build_sup_table_s9/…pol_proofreading.R", "Table S9", { outKind: "stale", note: "17 dominant signatures change", noteColor: HEX.accent2 });
  R("89 assignment, sample_info, Table S9", "build_sup_table_s10/", "Table S10", { outKind: "stale", note: "counts change", noteColor: HEX.accent2 });
  R("89 assignment, connection_table", "build_sup_table_s11/", "Table S11", { outKind: "stale", note: "819 tumors change", noteColor: HEX.accent2 });
  R("topography r.05, r.06, r.09 outputs", "build_sup_table_s12/ … s14/", "Tables S12–S14", { outKind: "stale", note: "needs new partial credit", noteColor: HEX.accent2 });
  R("build_fig_s6/*.txt", "build_sup_table_s15/, s16/", "Tables S15, S16", { inKind: "data", outKind: "stale", note: "hits change", noteColor: HEX.accent2 });
  R("sample_info, 83 spectra", "build_sup_table_s17/", "Table S17", { note: "unchanged" });
  legend(s, ["script", "data", "other", "output", "stale"]);
}

// ---- Slide 6: vignette -------------------------------------------------------
{
  pres.addSection({ title: "indelsig.net" });
  const s = newSlide(
    "Vignette: the indelsig.net signature pages",
    "vignette/ computes everything once, caches it, and renders one self-contained page per signature. Xueming puts the pages on the Shiny server.",
    "indelsig.net");
  box(s, "data", 0.5, 2.30, 3.2, 1.3, "signatures 83 / 89 / 476\nspectra 83 / 89 / 476\nconnection_table.tsv\nID476_ID89_mapping.txt\nper_sig_txt/", { fontSize: 10 });
  box(s, "script", 4.3, 2.60, 2.6, 0.7, "vignette.qmd", {});
  box(s, "data", 4.3, 3.80, 2.6, 0.7, "figure/standalone_data/\nall_sig_data.rds", { fontSize: 10 });
  box(s, "script", 7.5, 2.60, 2.6, 0.7, "render_separate_pages.R\n(onesig_standalone.qmd)", { fontSize: 10 });
  box(s, "script", 7.5, 3.80, 2.6, 0.7, "render_overview_table.R", { fontSize: 10 });
  box(s, "output", 10.6, 2.60, 2.25, 0.7, "44 signature pages", { bold: true });
  box(s, "output", 10.6, 3.80, 2.25, 0.7, "overview_table.html", { bold: true });
  box(s, "output", 4.3, 1.70, 2.6, 0.55, "combined HTML + PDF", { bold: true, fontSize: 10 });
  arrow(s, 3.7, 2.95, 4.3, 2.95);
  arrow(s, 5.6, 3.30, 5.6, 3.80);
  arrow(s, 5.6, 2.60, 5.6, 2.25);
  arrow(s, 6.9, 4.15, 7.5, 2.95); arrow(s, 6.9, 4.15, 7.5, 4.15);
  arrow(s, 10.1, 2.95, 10.6, 2.95); arrow(s, 10.1, 4.15, 10.6, 4.15);
  box(s, "output", 10.6, 4.90, 2.25, 0.6, "Shiny server copy\n(updated 2026-10-19)", { fontSize: 10 });
  arrow(s, 11.72, 4.50, 11.72, 4.90);
  s.addText([
    { text: "The vignette does not read the assignments. Rerendered with the new set on 2026-10-22: Table 1 and the pages have the same content.", options: { breakLine: true } },
    { text: "Rebuild: quarto render vignette/vignette.qmd, then the two render scripts, then send separate_pages/ to Xueming." },
  ], { x: 0.5, y: 5.75, w: 12.3, h: 0.7, fontSize: 12, color: HEX.dk2, isTextBox: true, margin: 0, paraSpaceAfter: 4 });
  legend(s, ["script", "data", "output"]);
}

// ---- Slide 7: what the new assignments change -------------------------------
{
  pres.addSection({ title: "New assignments" });
  const s = newSlide(
    "What Mo's assignments of 2026-10-19 change",
    "Rerun on 2026-10-22 in a separate worktree. The canonical assignment files and the paper are unchanged until we adopt the new set.",
    "New assignments");
  const stats = [
    ["819", "tumors with changed\nexposures (Mo: 784)"],
    ["731", "tumors fit with mSigAct\n(was 1,535)"],
    ["0", "Table S11 columns still\nsingle-counted (was 501)"],
  ];
  stats.forEach(([n, lab], i) => {
    const x = 0.5 + i * 2.2;
    s.addText(n, { x, y: 1.6, w: 2.0, h: 0.8, fontSize: 40, bold: true, color: HEX.accent3, fontFace: "Cambria", isTextBox: true, margin: 0 });
    s.addText(lab, { x, y: 2.4, w: 2.0, h: 0.6, fontSize: 11, color: HEX.dk2, isTextBox: true, margin: 0 });
  });
  const rows = [
    ["Figure 4A", "InsDel4 and InsDel10 in CNS and InsDel_N_beta in Lymphoid drop out; InsDel2c in Lymphoid appears"],
    ["Figure 4B", "Clusters reorder. Quoted correlations hardly move (InsDel1b–SBS10a 0.29 → 0.28)"],
    ["Figure S1", "Same 6,516 well-fit tumors as Mo's list; exposure color scale shifts"],
    ["Table S9", "Dominant signature changes for 17 of 248 tumors (InsDel7: 128 → 137)"],
    ["Table S10", "MSS InsDel_N_beta 2,629 → 2,409; InsDel1a without N_beta 2,596 → 2,724"],
    ["Table S15", "10 sex hits stay 10, but 5 instead of 6 pairs; Pancreas C_ID12/InsDel12 (q = 0.014) new; Other and Lymphoid hits gone"],
    ["Table S16", "Most age hits hold; InsDel4/C_ID4 now significant in breast and colon"],
    ["Topography", "Mo will redo partial credit; then Fig 5, S3, S4, Tables S12–S14"],
  ];
  let y = 1.55;
  for (const [what, change] of rows) {
    box(s, "stale", 7.0, y, 1.35, 0.5, what, { fontSize: 10, bold: true });
    s.addText(change, { x: 8.5, y, w: 4.35, h: 0.5, fontSize: 10, color: HEX.dk1, valign: "middle", isTextBox: true, margin: 0 });
    y += 0.6;
  }
  s.addText("Table S17 does not change. The methods text (5,440 + 1,535 tumors, routing rule) and the sex and age paragraphs must be revised when we adopt the new set.",
    { x: 0.5, y: 3.4, w: 6.1, h: 1.2, fontSize: 13, color: HEX.dk2, isTextBox: true, margin: 0 });
  legend(s, ["stale"]);
}

// ---- Slide 8: keeping ms.qmd in sync ---------------------------------------
{
  const s = newSlide(
    "Keeping the manuscript in sync with the code",
    "Numbers come from the built tables through ms_numbers.R, and three check scripts flag stale copies at render time.",
    "New assignments");
  box(s, "output", 0.5, 1.8, 3.0, 2.2, "Table S10\nTable S15\nfig4b_sbs_id_correlations.tsv\nbuild_fig5/*.cancer_proportion.xlsx\nsignatures, connection_table", { fontSize: 10 });
  box(s, "script", 4.1, 2.55, 2.6, 0.7, "ms_numbers.R\n(list ms_num)", {});
  box(s, "output", 7.3, 2.55, 2.4, 0.7, "inline `r ms_num$…`\nin ms.qmd", { fontSize: 11 });
  box(s, "output", 10.3, 2.55, 2.55, 0.7, "ms.pdf", { bold: true });
  arrow(s, 3.5, 2.9, 4.1, 2.9); arrow(s, 6.7, 2.9, 7.3, 2.9); arrow(s, 9.7, 2.9, 10.3, 2.9);

  box(s, "output", 0.5, 4.5, 3.0, 0.9, "built files in\nLiu2026_code_and_data", { fontSize: 11 });
  box(s, "script", 4.1, 4.5, 2.6, 0.9, "check_main_figures.R\ncheck_sup_figures.R\ncheck_sup_tables.R", { fontSize: 10 });
  box(s, "other", 4.1, 5.6, 2.6, 0.55, "*_manifest.tsv (paper repo)", { fontSize: 10 });
  box(s, "output", 10.3, 4.5, 2.55, 0.9, "warnings in ms.pdf\nwhen a copy is stale", { fontSize: 10 });
  arrow(s, 3.5, 4.95, 4.1, 4.95); arrow(s, 5.4, 5.6, 5.4, 5.4); arrow(s, 6.7, 4.95, 10.3, 4.95);
  s.addText("Set LIU2026_CODE_REPO to read the tables from another checkout, for example the worktree built with the new assignments.",
    { x: 7.3, y: 5.55, w: 5.55, h: 0.7, fontSize: 12, color: HEX.dk2, isTextBox: true, margin: 0 });
  legend(s, ["script", "other", "output"]);
}

(async () => {
  const out = process.argv[2] || "data_and_analysis_flow.pptx";
  await pres.writeFile({ fileName: out });
  await applyTheme(out, THEME);
  console.log("wrote", out);
})();
