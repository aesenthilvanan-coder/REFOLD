import Link from "next/link";
import { PCDEntry } from "../../types";
import { notFound } from "next/navigation";
import { MolViewer } from "../../components/MolViewer";
import { fetchAtlas } from "../../lib/atlas";

export const dynamic = "force-dynamic";

function ScoreRow({ label, value, color = "var(--accent)", max = 1 }: {
  label: string; value: number; color?: string; max?: number;
}) {
  const pct = Math.min((value / max) * 100, 100);
  return (
    <div>
      <div className="flex justify-between items-baseline mb-1.5">
        <span className="section-label">{label}</span>
        <span className="text-xs font-mono font-bold" style={{ color }}>{value.toFixed(3)}</span>
      </div>
      <div className="score-bar-track">
        <div className="score-bar-fill" style={{ width: `${pct}%`, background: color }} />
      </div>
    </div>
  );
}

function PropCell({ label, value, pass }: { label: string; value: string; pass?: boolean }) {
  const color = pass === true ? "var(--green)" : pass === false ? "var(--red)" : "var(--accent)";
  return (
    <div className="card p-3 text-center space-y-1">
      <div className="text-sm font-bold font-mono" style={{ color }}>{value}</div>
      <div className="data-label">{label}</div>
    </div>
  );
}

function EijHeatmap({ labels, values }: { labels: string[]; values: number[][] }) {
  if (!values?.length) return (
    <div className="card p-6 text-center">
      <span className="section-label">Matrix data unavailable</span>
    </div>
  );

  const maxVal = Math.max(...values.flat().filter(v => v > 0));

  function cellBg(v: number) {
    if (v === 0) return "var(--bg)";
    const t = Math.min(v / maxVal, 1);
    const r = Math.round(6 + t * 133);
    const g = Math.round(182 - t * 90);
    const b = Math.round(212 + t * 34);
    return `rgba(${r},${g},${b},0.7)`;
  }

  return (
    <div className="overflow-x-auto">
      <table className="text-[10px] font-mono border-collapse">
        <thead>
          <tr>
            <th className="p-1" />
            {labels.map(l => (
              <th key={l} className="p-1 text-center font-semibold" style={{ color: "var(--text-secondary)", minWidth: 38 }}>{l}</th>
            ))}
          </tr>
        </thead>
        <tbody>
          {values.map((row, i) => (
            <tr key={i}>
              <td className="pr-2 font-semibold text-right whitespace-nowrap" style={{ color: "var(--text-secondary)" }}>{labels[i]}</td>
              {row.map((v, j) => (
                <td key={j} className="p-0.5">
                  <div className="flex items-center justify-center rounded"
                       style={{ background: cellBg(v), width: 38, height: 26, opacity: i === j ? 0.2 : 1 }}
                       title={`${labels[i]}↔${labels[j]}: ${v.toFixed(2)} Å`}>
                    <span style={{ color: "white", fontSize: 9, fontWeight: 700 }}>{v.toFixed(1)}</span>
                  </div>
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
      <div className="flex items-center gap-2 mt-3">
        <span className="text-[10px]" style={{ color: "var(--text-muted)" }}>0 Å</span>
        <div className="flex-1 h-1.5 rounded-full" style={{ background: "linear-gradient(90deg, var(--bg), var(--accent), #8b5cf6)" }} />
        <span className="text-[10px]" style={{ color: "var(--text-muted)" }}>{maxVal.toFixed(0)} Å</span>
      </div>
    </div>
  );
}

export default async function EntryPage({ params }: { params: Promise<{ id: string }> }) {
  const { id } = await params;
  const atlas = await fetchAtlas();
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const entry: any = atlas.entries.find((e: any) => e.entry_id === id);
  if (!entry) notFound();

  const pocket   = entry.pocket   ?? {};
  const chap     = entry.chaperone ?? {};
  const seq      = entry.sequence  ?? {};
  const eij      = entry.eij_matrix ?? {};
  const assets   = entry.assets    ?? {};
  const meta     = entry.metadata  ?? {};

  const pocketLabelList: string[] = (seq.pocket_lining_residues ?? "").split("-").filter(Boolean);
  const pdbUrl = entry.pdb_structure ?? assets.pdb_structure ?? "";
  // binding_mode: two possible formats from pipeline
  //   new entries:  { salt_bridge: "description...", ... }  (key=type, value=description)
  //   backfilled:   { type: "salt_bridge", residues: [...] }
  const bindingModeRaw = chap.binding_mode ?? {};
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bindingModeEntries: [string, any][] = "type" in bindingModeRaw
    ? [[bindingModeRaw.type ?? "interaction", Array.isArray(bindingModeRaw.residues) ? bindingModeRaw.residues.join(", ") : ""]]
    : Object.entries(bindingModeRaw);

  return (
    <div style={{ background: "var(--bg)", minHeight: "100vh" }}>
      {/* Nav */}
      <nav className="fixed top-0 left-0 right-0 z-50 border-b border-[var(--border)] backdrop-blur-md"
           style={{ background: "rgba(7,9,15,0.92)" }}>
        <div className="max-w-7xl mx-auto px-6 flex items-center gap-2 h-14 text-sm">
          <Link href="/" className="flex items-center gap-2">
            <div className="w-6 h-6 rounded border border-[var(--accent-border)] flex items-center justify-center"
                 style={{ background: "var(--accent-dim)" }}>
              <span className="font-bold text-[10px] font-mono" style={{ color: "var(--accent)" }}>P</span>
            </div>
            <span style={{ color: "var(--text-secondary)" }}>PCD</span>
          </Link>
          <span style={{ color: "var(--text-muted)" }}>/</span>
          <Link href="/database" style={{ color: "var(--text-secondary)" }} className="hover:text-white transition-colors">Database</Link>
          <span style={{ color: "var(--text-muted)" }}>/</span>
          <span className="font-mono font-semibold text-xs" style={{ color: "var(--text-primary)" }}>{entry.entry_id}</span>
        </div>
      </nav>

      <main className="max-w-7xl mx-auto px-6 pt-20 pb-16 space-y-6">

        {/* ── Header ─────────────────────────────────────────────────── */}
        <div className="card rounded-xl p-6 grid lg:grid-cols-3 gap-6">
          <div className="lg:col-span-2 space-y-4">
            <div className="flex flex-wrap items-center gap-2">
              <span className="tag tag-cyan">{meta.gene}</span>
              <span className="tag tag-amber">{meta.mutation_mature}</span>
              <span className="tag tag-violet">{meta.variant_class ?? "Pathogenic"}</span>
              <span className="tag tag-green">{entry.status ?? "COMPLETE"}</span>
            </div>
            <h1 className="text-2xl font-black" style={{ color: "var(--text-primary)" }}>{meta.disease}</h1>
            <p className="text-sm leading-6" style={{ color: "var(--text-secondary)" }}>{meta.mechanism ?? "Protein misfolding pathogenic variant"}</p>
            <div className="flex flex-wrap gap-4 text-xs font-mono" style={{ color: "var(--text-muted)" }}>
              <span>UniProt <span style={{ color: "var(--text-secondary)" }}>{meta.uniprot}</span></span>
              <span>ClinVar <span style={{ color: "var(--text-secondary)" }}>#{meta.clinvar_id}</span></span>
              <span>Precursor <span style={{ color: "var(--text-secondary)" }}>{meta.mutation_precursor ?? meta.mutation_mature}</span></span>
            </div>
          </div>
          <div className="space-y-4">
            <div className="section-label">Key Scores</div>
            <ScoreRow label="fpocket Druggability" value={pocket.fpocket_druggability ?? 0} color="var(--accent)" />
            <ScoreRow label="Chaperone Composite" value={chap.composite_score ?? 0} color="var(--violet)" />
            <ScoreRow label="Pocket Affinity" value={chap.pocket_affinity_score ?? chap.composite_score ?? 0} color="#3b82f6" />
            <ScoreRow label="Drug-likeness (QED)" value={chap.qed ?? 0} color="var(--green)" />
          </div>
        </div>

        {/* ── Pocket stats + 3D viewer side by side ──────────────────── */}
        <div className="card rounded-xl p-5">
          <div className="grid lg:grid-cols-[1fr_340px] gap-6">

            {/* Left: pocket stats */}
            <div className="space-y-4">
              <div>
                <span className="tag tag-cyan">Stage 2 — Transient Pocket</span>
                <h2 className="mt-3 font-bold" style={{ color: "var(--text-primary)" }}>Cryptic Binding Site</h2>
              </div>

              <div className="grid grid-cols-2 gap-2">
                {[
                  { l: "fpocket Drug.", v: (pocket.fpocket_druggability ?? 0).toFixed(3), c: "var(--accent)" },
                  { l: "WT Baseline", v: (pocket.wt_baseline_druggability ?? 0).toFixed(3), c: "var(--red)" },
                  { l: "Volume (Å³)", v: (pocket.volume_angstrom3 ?? 0).toFixed(1), c: "var(--text-secondary)" },
                  { l: "α-Spheres", v: String(pocket.alpha_sphere_count ?? "—"), c: "var(--text-secondary)" },
                  { l: "Conformations", v: String(pocket.n_conformations_sampled ?? 20), c: "var(--text-secondary)" },
                  { l: "Mut→Pocket Dist.", v: pocket.dist_mutation_to_pocket_angstrom != null ? `${(pocket.dist_mutation_to_pocket_angstrom).toFixed(1)} Å` : "—", c: "var(--amber)" },
                ].map(p => (
                  <div key={p.l} className="rounded-lg border border-[var(--border)] p-3"
                       style={{ background: "var(--bg-1)" }}>
                    <div className="text-sm font-bold font-mono" style={{ color: p.c }}>{p.v}</div>
                    <div className="data-label mt-0.5">{p.l}</div>
                  </div>
                ))}
              </div>

              <div className="grid grid-cols-2 gap-2">
                <div className="rounded-lg border border-[var(--border)] p-3" style={{ background: "var(--bg-1)" }}>
                  <div className="data-label mb-1">Pocket Type</div>
                  <div className="text-xs font-semibold" style={{ color: "var(--text-secondary)" }}>
                    {pocket.pocket_type ?? (pocket.cryptic ? "Cryptic" : "Allosteric")}
                  </div>
                </div>
                <div className="rounded-lg border border-[var(--border)] p-3" style={{ background: "var(--bg-1)" }}>
                  <div className="data-label mb-1">Pocket Center (Å)</div>
                  <div className="font-mono text-[11px]" style={{ color: "var(--accent)" }}>
                    ({(pocket.center_angstrom ?? [0,0,0]).map((v: number) => v.toFixed(1)).join(", ")})
                  </div>
                </div>
              </div>

              <div>
                <div className="data-label mb-2">Pocket-Lining Residues</div>
                <div className="flex flex-wrap gap-1.5">
                  {pocketLabelList.map(r => (
                    <span key={r} className="tag tag-violet text-[10px]">{r}</span>
                  ))}
                </div>
              </div>

              <div className="rounded-lg border border-[var(--border)] p-3" style={{ background: "var(--bg-1)" }}>
                <div className="data-label mb-1">Sequence Slice</div>
                <div className="font-mono text-[11px] break-all leading-5" style={{ color: "var(--text-secondary)" }}>
                  {seq.sequence_slice_around_pocket ?? "—"}
                </div>
              </div>
            </div>

            {/* Right: square 3D viewer */}
            <div className="space-y-2">
              <div className="flex items-center gap-2">
                <div className="section-label">3D Structure</div>
                <div className="w-1.5 h-1.5 rounded-full bg-[var(--accent)]" />
                <span className="text-[10px]" style={{ color: "var(--text-muted)" }}>drag to rotate</span>
              </div>
              {pdbUrl ? (
                <MolViewer
                  pdbUrl={pdbUrl}
                  pocketResidues={pocketLabelList}
                  height={340}
                />
              ) : (
                <div className="rounded-xl border border-[var(--border)] flex items-center justify-center"
                     style={{ height: 340, background: "var(--bg-1)" }}>
                  <span className="section-label">No structure</span>
                </div>
              )}
              <p className="text-[10px] leading-4" style={{ color: "var(--text-muted)" }}>
                {pocket.target_conformation ?? "ANM-sampled transient conformation"}
              </p>
            </div>
          </div>
        </div>

        {/* ── Chaperone ──────────────────────────────────────────────── */}
        <div className="card rounded-xl p-5 space-y-5">
            <div>
              <span className="tag tag-violet">Stage 3 — De Novo Chaperone</span>
              <h2 className="mt-3 font-bold" style={{ color: "var(--text-primary)" }}>{chap.common_name ?? "De Novo Candidate"}</h2>
              <p className="text-xs italic mt-1" style={{ color: "var(--text-secondary)" }}>{chap.iupac_name ?? chap.smiles ?? "—"}</p>
            </div>

            <div className="rounded-lg border border-[var(--border)] p-3" style={{ background: "var(--bg-1)" }}>
              <div className="data-label mb-2">SMILES</div>
              <div className="font-mono text-xs break-all" style={{ color: "var(--accent)" }}>
                {chap.smiles ?? "—"}
              </div>
            </div>

            <div>
              <div className="data-label mb-2">Physicochemical Properties</div>
              <div className="grid grid-cols-4 gap-1.5">
                <PropCell label="MW (Da)" value={(chap.mw ?? 0).toFixed(0)} />
                <PropCell label="logP" value={(chap.logp ?? 0).toFixed(2)} pass={(chap.logp ?? 0) >= 0 && (chap.logp ?? 0) <= 5} />
                <PropCell label="QED" value={(chap.qed ?? 0).toFixed(3)} pass={(chap.qed ?? 0) > 0.6} />
                <PropCell label="SA" value={(chap.sa_score ?? 0).toFixed(1)} pass={(chap.sa_score ?? 0) < 4} />
                <PropCell label="HBD" value={String(chap.hbd ?? 0)} pass={(chap.hbd ?? 0) <= 5} />
                <PropCell label="HBA" value={String(chap.hba ?? 0)} pass={(chap.hba ?? 0) <= 10} />
                <PropCell label="TPSA" value={(chap.tpsa ?? 0).toFixed(0)} pass={(chap.tpsa ?? 0) <= 140} />
              </div>
            </div>

            <div>
              <div className="data-label mb-2">Drug-likeness</div>
              <div className="flex flex-wrap gap-2">
                {[
                  { label: "Lipinski Ro5", pass: typeof chap.lipinski === "object" ? chap.lipinski?.passes : (chap.lipinski ?? ((chap.mw ?? 0) <= 500 && (chap.hbd ?? 0) <= 5 && (chap.hba ?? 0) <= 10 && (chap.logp ?? 0) <= 5)) },
                  { label: "Veber Oral", pass: typeof chap.veber === "object" ? chap.veber?.passes : (chap.veber ?? ((chap.rotatable_bonds ?? 0) <= 10 && (chap.tpsa ?? 0) <= 140)) },
                  { label: "PAINS-free", pass: chap.pains != null ? !chap.pains : true },
                ].map(r => (
                  <div key={r.label} className="flex items-center gap-1.5 px-2.5 py-1.5 rounded border text-xs font-semibold"
                       style={{
                         background: r.pass ? "var(--green-dim)" : "var(--red-dim)",
                         borderColor: r.pass ? "#10b98130" : "#ef444430",
                         color: r.pass ? "var(--green)" : "var(--red)",
                       }}>
                    {r.pass ? "✓" : "✗"} {r.label}
                  </div>
                ))}
              </div>
            </div>

            <div>
              <div className="data-label mb-2">Predicted Binding Mode</div>
              <div className="space-y-1.5">
                {bindingModeEntries.map(([k, v]) => (
                  <div key={k} className="flex gap-2.5 rounded-lg border border-[var(--border)] p-2.5"
                       style={{ background: "var(--bg-1)" }}>
                    <div className="w-1 h-1 rounded-full mt-1.5 shrink-0" style={{ background: "var(--violet)" }} />
                    <div>
                      <div className="text-[10px] font-semibold uppercase tracking-wider mb-0.5"
                           style={{ color: "var(--violet)" }}>
                        {String(k).replace(/_/g, " ")}
                      </div>
                      <div className="text-xs leading-5" style={{ color: "var(--text-secondary)" }}>
                        {Array.isArray(v) ? v.join(", ") : String(v)}
                      </div>
                    </div>
                  </div>
                ))}
              </div>
            </div>
        </div>

        {/* ── E_ij Matrix ─────────────────────────────────────────────── */}
        <div className="card rounded-xl p-5 space-y-4">
          <div>
            <span className="tag tag-amber">Stage 2 Output</span>
            <h2 className="mt-3 font-bold" style={{ color: "var(--text-primary)" }}>
              E<sub>ij</sub> Pairwise Distance Matrix
            </h2>
            <p className="text-sm mt-1 leading-6" style={{ color: "var(--text-secondary)" }}>
              Cα–Cα distances (Å) between pocket-lining residues in the selected transient conformation.
              Used as conditioning input for Stage 3 pharmacophore-guided SMILES generation.
            </p>
          </div>
          <div className="flex gap-4 text-xs font-mono" style={{ color: "var(--text-muted)" }}>
            <span>Shape <span style={{ color: "var(--text-secondary)" }}>{(eij.shape ?? [0,0])[0]} × {(eij.shape ?? [0,0])[1]}</span></span>
            <span>Source <span style={{ color: "var(--text-secondary)" }}>{eij.file ?? "—"}</span></span>
          </div>
          <EijHeatmap labels={eij.residue_labels ?? []} values={eij.values ?? []} />
        </div>

        {/* ── Sequence ─────────────────────────────────────────────────── */}
        <div className="card rounded-xl p-5 space-y-4">
          <div>
            <span className="tag tag-green">Sequence</span>
            <p className="mt-3 text-xs font-mono" style={{ color: "var(--text-muted)" }}>{seq.fasta_header ?? `>${meta.gene ?? "?"} ${meta.mutation_mature ?? ""}`}</p>
          </div>
          <div className="rounded-lg border border-[var(--border)] p-4 overflow-x-auto"
               style={{ background: "var(--bg-1)" }}>
            <div className="font-mono text-[11px] break-all leading-5">
              {(seq.fasta_sequence ?? seq.full_sequence ?? "").split("").map((ch: string, i: number) => {
                const pos = i + 1;
                const isLining = seq.pocket_lining_positions_precursor?.includes(pos);
                return (
                  <span key={i}
                        style={{ color: isLining ? "var(--accent)" : "var(--text-muted)", fontWeight: isLining ? 700 : 400 }}
                        title={isLining ? `Pos ${pos}` : undefined}>
                    {ch}
                  </span>
                );
              })}
            </div>
            <div className="mt-3 flex items-center gap-2">
              <span className="inline-block w-3 h-2 rounded-sm" style={{ background: "var(--accent)" }} />
              <span className="text-[10px]" style={{ color: "var(--text-muted)" }}>Pocket-lining positions</span>
            </div>
          </div>
          <div className="text-xs font-mono" style={{ color: "var(--text-muted)" }}>
            1-letter: <span className="font-bold tracking-widest" style={{ color: "var(--text-secondary)" }}>
              {seq.pocket_lining_1letter ?? "—"}
            </span>
          </div>
        </div>

        {/* ── Assets ──────────────────────────────────────────────────── */}
        <div className="card rounded-xl p-5 space-y-4">
          <span className="tag tag-cyan">Data Assets</span>
          <div className="grid sm:grid-cols-2 gap-2">
            {[
              { label: "E_ij Matrix (.npy)", path: assets.eij_matrix_npy },
              ...(assets.eij_matrix_csv ? [{ label: "E_ij Matrix (.csv)", path: assets.eij_matrix_csv }] : []),
              { label: "Transient Conformation (.pdb)", path: assets.transient_conformation_pdb },
              { label: "Pocket Summary (.json)", path: assets.pocket_summary_json },
              { label: "Stage 3 Candidates (.json)", path: assets.stage3_candidates_json },
            ].map(a => (
              <div key={a.label} className="rounded-lg border border-[var(--border)] p-3"
                   style={{ background: "var(--bg-1)" }}>
                <div className="data-label mb-1">{a.label}</div>
                <div className="font-mono text-[11px] break-all" style={{ color: "var(--text-secondary)" }}>{a.path}</div>
              </div>
            ))}
          </div>
        </div>

        {/* Footer nav */}
        <div className="flex justify-between items-center pt-2">
          <Link href="/database"
                className="px-4 py-2 rounded-lg text-sm border transition-colors"
                style={{ borderColor: "var(--border)", color: "var(--text-secondary)", background: "var(--bg-2)" }}>
            ← Back to Database
          </Link>
          <div className="text-[10px] font-mono" style={{ color: "var(--text-muted)" }}>
            REFOLD · S.Y.A.L.I.S Labs · {entry.investigator}
          </div>
        </div>

      </main>
    </div>
  );
}
