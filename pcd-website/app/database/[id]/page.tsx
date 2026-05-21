import Link from "next/link";
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
      <div className="flex justify-between items-baseline mb-1">
        <span className="data-label">{label}</span>
        <span className="font-mono font-bold tabular-nums" style={{ color, fontSize: 12 }}>{value.toFixed(3)}</span>
      </div>
      <div className="score-bar-track">
        <div className="score-bar-fill" style={{ width: `${pct}%`, background: color }} />
      </div>
    </div>
  );
}

function PropCell({ label, value, pass }: { label: string; value: string; pass?: boolean }) {
  const color = pass === true
    ? "var(--green)"
    : pass === false
    ? "var(--red)"
    : "var(--text-primary)";
  const bg = pass === true
    ? "var(--green-dim)"
    : pass === false
    ? "var(--red-dim)"
    : "var(--bg-1)";
  return (
    <div className="p-2.5 text-center border" style={{ background: bg, borderColor: "var(--border)", borderRadius: "2px" }}>
      <div className="text-sm font-bold font-mono tabular-nums" style={{ color }}>{value}</div>
      <div className="data-label mt-0.5">{label}</div>
    </div>
  );
}

function EijHeatmap({ labels, values }: { labels: string[]; values: number[][] }) {
  if (!values?.length) return (
    <div className="p-6 text-center border" style={{ borderColor: "var(--border)", borderRadius: "3px" }}>
      <span className="section-label">Matrix data unavailable</span>
    </div>
  );

  const maxVal = Math.max(...values.flat().filter(v => v > 0));

  function cellBg(v: number) {
    if (v === 0) return "var(--bg-2)";
    const t = Math.min(v / maxVal, 1);
    // White → light blue → deep navy
    const opacity = 0.08 + t * 0.88;
    return `rgba(30, 58, 138, ${opacity.toFixed(2)})`;
  }

  function cellText(v: number) {
    if (v === 0) return "var(--text-muted)";
    const t = Math.min(v / maxVal, 1);
    return t > 0.42 ? "#FFFFFF" : "var(--text-secondary)";
  }

  return (
    <div className="overflow-x-auto">
      <table className="text-[10px] font-mono border-collapse">
        <thead>
          <tr>
            <th className="p-1" />
            {labels.map(l => (
              <th key={l} className="p-1 text-center font-semibold"
                  style={{ color: "var(--text-muted)", minWidth: 38 }}>{l}</th>
            ))}
          </tr>
        </thead>
        <tbody>
          {values.map((row, i) => (
            <tr key={i}>
              <td className="pr-2 font-semibold text-right whitespace-nowrap"
                  style={{ color: "var(--text-muted)" }}>{labels[i]}</td>
              {row.map((v, j) => (
                <td key={j} className="p-0.5">
                  <div className="flex items-center justify-center"
                       style={{
                         background: cellBg(v),
                         width: 38, height: 26,
                         opacity: i === j ? 0.15 : 1,
                       }}
                       title={`${labels[i]}↔${labels[j]}: ${v.toFixed(2)} Å`}>
                    <span style={{ color: cellText(v), fontSize: 9, fontWeight: 700 }}>
                      {v.toFixed(1)}
                    </span>
                  </div>
                </td>
              ))}
            </tr>
          ))}
        </tbody>
      </table>
      {/* Legend */}
      <div className="flex items-center gap-2 mt-3">
        <span className="text-[10px]" style={{ color: "var(--text-muted)" }}>0 Å</span>
        <div className="flex-1 h-2 flex border" style={{ borderColor: "var(--border)" }}>
          {[0.08,0.25,0.42,0.59,0.76,0.96].map(t => (
            <div key={t} className="flex-1" style={{ background: `rgba(30,58,138,${t})` }} />
          ))}
        </div>
        <span className="text-[10px]" style={{ color: "var(--text-muted)" }}>{maxVal.toFixed(0)} Å</span>
      </div>
    </div>
  );
}

// ── Section wrapper ──────────────────────────────────────────────────────────
function Panel({ children, className = "" }: { children: React.ReactNode; className?: string }) {
  return (
    <div className={`border p-5 ${className}`}
         style={{ borderColor: "var(--border)", borderRadius: "3px", background: "var(--bg)" }}>
      {children}
    </div>
  );
}

function PanelHeading({ tag, tagClass = "tag-cyan", title }: { tag: string; tagClass?: string; title: string }) {
  return (
    <div className="mb-4">
      <span className={`tag ${tagClass}`}>{tag}</span>
      <h2 className="mt-2 text-base font-bold" style={{ color: "var(--text-primary)" }}>{title}</h2>
    </div>
  );
}

// ── Stat cell (2×3 grid) ─────────────────────────────────────────────────────
function StatCell({ label, value, color = "var(--text-primary)" }: { label: string; value: string; color?: string }) {
  return (
    <div className="p-3 border" style={{ background: "var(--bg-1)", borderColor: "var(--border)", borderRadius: "2px" }}>
      <div className="font-mono font-bold tabular-nums text-sm" style={{ color }}>{value}</div>
      <div className="data-label mt-0.5">{label}</div>
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

  const bindingModeRaw = chap.binding_mode ?? {};
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bindingModeEntries: [string, any][] = "type" in bindingModeRaw
    ? [[bindingModeRaw.type ?? "interaction", Array.isArray(bindingModeRaw.residues) ? bindingModeRaw.residues.join(", ") : ""]]
    : Object.entries(bindingModeRaw);

  return (
    <div style={{ background: "var(--bg-1)", minHeight: "100vh" }}>

      {/* ── Nav ── */}
      <nav className="fixed top-0 left-0 right-0 z-50 border-b"
           style={{ background: "var(--bg)", borderColor: "var(--border)" }}>
        <div className="max-w-7xl mx-auto px-6 flex items-center gap-2 h-10 text-[11px]">
          <Link href="/" className="font-mono font-bold tracking-widest"
                style={{ color: "var(--accent)" }}>PCD</Link>
          <span style={{ color: "var(--border)" }}>|</span>
          <Link href="/database" className="nav-link">Database</Link>
          <span style={{ color: "var(--border-mid)" }}>/</span>
          <span className="font-mono font-semibold text-[11px]" style={{ color: "var(--text-primary)" }}>{entry.entry_id}</span>
        </div>
      </nav>

      <main className="max-w-7xl mx-auto px-6 pb-12 space-y-4" style={{ paddingTop: "52px" }}>

        {/* ── Header ── */}
        <Panel className="!p-0 overflow-hidden">
          <div className="px-5 py-4 border-b" style={{ borderColor: "var(--border)", background: "var(--bg-1)" }}>
            <div className="flex flex-wrap items-center gap-2">
              <span className="tag tag-cyan">{meta.gene}</span>
              <span className="tag tag-amber">{meta.mutation_mature}</span>
              <span className="tag tag-violet">{meta.variant_class ?? "Pathogenic"}</span>
              <span className="tag tag-green">{entry.status ?? "COMPLETE"}</span>
              <span className="font-mono text-[10px] ml-2" style={{ color: "var(--text-muted)" }}>{entry.entry_id}</span>
            </div>
            <h1 className="mt-2 text-xl font-bold" style={{ color: "var(--text-primary)" }}>{meta.disease}</h1>
            <p className="text-xs mt-1 leading-5" style={{ color: "var(--text-secondary)" }}>{meta.mechanism ?? "Protein misfolding pathogenic variant"}</p>
            <div className="flex flex-wrap gap-4 mt-2 text-[11px] font-mono" style={{ color: "var(--text-muted)" }}>
              <span>UniProt <span style={{ color: "var(--text-secondary)" }}>{meta.uniprot}</span></span>
              <span>ClinVar <span style={{ color: "var(--text-secondary)" }}>#{meta.clinvar_id}</span></span>
              <span>Precursor <span style={{ color: "var(--text-secondary)" }}>{meta.mutation_precursor ?? meta.mutation_mature}</span></span>
            </div>
          </div>
          <div className="px-5 py-4 grid sm:grid-cols-4 gap-4">
            <ScoreRow label="fpocket Druggability" value={pocket.fpocket_druggability ?? 0} color="var(--accent)" />
            <ScoreRow label="Chaperone Composite" value={chap.composite_score ?? 0} color="var(--violet)" />
            <ScoreRow label="Pocket Affinity" value={chap.pocket_affinity_score ?? chap.composite_score ?? 0} color="#2563EB" />
            <ScoreRow label="Drug-likeness (QED)" value={chap.qed ?? 0} color="var(--green)" />
          </div>
        </Panel>

        {/* ── Pocket + 3D viewer ── */}
        <div className="grid lg:grid-cols-[1fr_320px] gap-4">

          <Panel>
            <PanelHeading tag="Stage 2 — Transient Pocket" title="Cryptic Binding Site" />

            {/* Stats grid */}
            <div className="grid grid-cols-2 gap-2 mb-4">
              <StatCell label="fpocket Drug." value={(pocket.fpocket_druggability ?? 0).toFixed(3)} color="var(--accent)" />
              <StatCell label="WT Baseline" value={(pocket.wt_baseline_druggability ?? 0).toFixed(3)} color="var(--red)" />
              <StatCell label="Volume (Å³)" value={(pocket.volume_angstrom3 ?? 0).toFixed(1)} />
              <StatCell label="α-Spheres" value={String(pocket.alpha_sphere_count ?? "—")} />
              <StatCell label="Conformations" value={String(pocket.n_conformations_sampled ?? 20)} />
              <StatCell label="Mut→Pocket Dist."
                value={pocket.dist_mutation_to_pocket_angstrom != null
                  ? `${pocket.dist_mutation_to_pocket_angstrom.toFixed(1)} Å` : "—"}
                color="var(--amber)" />
            </div>

            {/* Pocket type + center */}
            <div className="grid grid-cols-2 gap-2 mb-4">
              <div className="p-3 border" style={{ background: "var(--bg-1)", borderColor: "var(--border)", borderRadius: "2px" }}>
                <div className="data-label mb-1">Pocket Type</div>
                <div className="text-xs font-semibold" style={{ color: "var(--text-secondary)" }}>
                  {pocket.pocket_type ?? (pocket.cryptic ? "Cryptic" : "Allosteric")}
                </div>
              </div>
              <div className="p-3 border" style={{ background: "var(--bg-1)", borderColor: "var(--border)", borderRadius: "2px" }}>
                <div className="data-label mb-1">Pocket Center (Å)</div>
                <div className="font-mono text-[11px]" style={{ color: "var(--accent)" }}>
                  ({(pocket.center_angstrom ?? [0,0,0]).map((v: number) => v.toFixed(1)).join(", ")})
                </div>
              </div>
            </div>

            {/* Pocket-lining residues */}
            <div>
              <div className="data-label mb-2">Pocket-Lining Residues</div>
              <div className="flex flex-wrap gap-1">
                {pocketLabelList.map(r => (
                  <span key={r} className="tag tag-violet text-[9px]">{r}</span>
                ))}
              </div>
            </div>

            {/* Sequence slice */}
            <div className="mt-4 p-3 border" style={{ background: "var(--bg-1)", borderColor: "var(--border)", borderRadius: "2px" }}>
              <div className="data-label mb-1">Sequence Slice (±20 aa around pocket)</div>
              <div className="font-mono text-[11px] break-all leading-5" style={{ color: "var(--text-secondary)" }}>
                {seq.sequence_slice_around_pocket ?? "—"}
              </div>
            </div>
          </Panel>

          {/* 3D Viewer */}
          <Panel className="flex flex-col gap-3">
            <div>
              <div className="section-label mb-1">3D Structure</div>
              <div className="text-[10px]" style={{ color: "var(--text-muted)" }}>Transient pocket conformation · drag to rotate</div>
            </div>
            {pdbUrl ? (
              <MolViewer
                pdbUrl={pdbUrl}
                pocketResidues={pocketLabelList}
                ligandSdfUrl={assets.ligand_sdf_url ?? undefined}
                height={300}
              />
            ) : (
              <div className="flex items-center justify-center border"
                   style={{ height: 300, background: "var(--bg-1)", borderColor: "var(--border)", borderRadius: "2px" }}>
                <span className="section-label">Structure unavailable</span>
              </div>
            )}
            <p className="text-[10px] leading-4" style={{ color: "var(--text-muted)" }}>
              {pocket.target_conformation ?? "ANM-sampled transient conformation"}
            </p>
          </Panel>
        </div>

        {/* ── Chaperone ── */}
        <Panel>
          <PanelHeading tag="Stage 3 — De Novo Chaperone" tagClass="tag-violet"
                        title={chap.common_name ?? "De Novo Candidate"} />

          <div className="space-y-4">
            {/* SMILES */}
            <div className="p-3 border" style={{ background: "var(--bg-1)", borderColor: "var(--border)", borderRadius: "2px" }}>
              <div className="data-label mb-1">SMILES String</div>
              <div className="font-mono text-xs break-all leading-5" style={{ color: "var(--accent)" }}>
                {chap.smiles ?? "—"}
              </div>
            </div>
            {chap.iupac_name && (
              <div className="text-xs italic" style={{ color: "var(--text-muted)" }}>{chap.iupac_name}</div>
            )}

            {/* Properties grid */}
            <div>
              <div className="data-label mb-2">Physicochemical Properties</div>
              <div className="grid grid-cols-4 gap-1.5">
                <PropCell label="MW (Da)" value={(chap.mw ?? 0).toFixed(0)} />
                <PropCell label="logP" value={(chap.logp ?? 0).toFixed(2)}
                          pass={(chap.logp ?? 0) >= 0 && (chap.logp ?? 0) <= 5} />
                <PropCell label="QED" value={(chap.qed ?? 0).toFixed(3)}
                          pass={(chap.qed ?? 0) > 0.6} />
                <PropCell label="SA" value={(chap.sa_score ?? 0).toFixed(1)}
                          pass={(chap.sa_score ?? 0) < 4} />
                <PropCell label="HBD" value={String(chap.hbd ?? 0)} pass={(chap.hbd ?? 0) <= 5} />
                <PropCell label="HBA" value={String(chap.hba ?? 0)} pass={(chap.hba ?? 0) <= 10} />
                <PropCell label="TPSA (Å²)" value={(chap.tpsa ?? 0).toFixed(0)} pass={(chap.tpsa ?? 0) <= 140} />
              </div>
            </div>

            {/* Drug-likeness flags */}
            <div>
              <div className="data-label mb-2">Drug-likeness Rules</div>
              <div className="flex flex-wrap gap-2">
                {[
                  { label: "Lipinski Ro5", pass: typeof chap.lipinski === "object" ? chap.lipinski?.passes : (chap.lipinski ?? ((chap.mw ?? 0) <= 500 && (chap.hbd ?? 0) <= 5 && (chap.hba ?? 0) <= 10 && (chap.logp ?? 0) <= 5)) },
                  { label: "Veber Oral", pass: typeof chap.veber === "object" ? chap.veber?.passes : (chap.veber ?? ((chap.rotatable_bonds ?? 0) <= 10 && (chap.tpsa ?? 0) <= 140)) },
                  { label: "PAINS-free", pass: chap.pains != null ? !chap.pains : true },
                ].map(r => (
                  <div key={r.label}
                       className="flex items-center gap-1.5 px-3 py-1.5 border text-xs font-semibold"
                       style={{
                         background: r.pass ? "var(--green-dim)" : "var(--red-dim)",
                         borderColor: r.pass ? "var(--green-border)" : "var(--red-border)",
                         color: r.pass ? "var(--green)" : "var(--red)",
                         borderRadius: "2px",
                       }}>
                    {r.pass ? "✓" : "✗"} {r.label}
                  </div>
                ))}
              </div>
            </div>

            {/* Binding mode */}
            <div>
              <div className="data-label mb-2">Predicted Binding Mode</div>
              <div className="space-y-1.5">
                {bindingModeEntries.map(([k, v]) => (
                  <div key={k} className="flex gap-2.5 border p-2.5"
                       style={{ background: "var(--bg-1)", borderColor: "var(--border)", borderRadius: "2px" }}>
                    <div className="w-0.5 shrink-0 rounded-sm mt-0.5" style={{ background: "var(--violet)" }} />
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
        </Panel>

        {/* ── E_ij Matrix ── */}
        <Panel>
          <PanelHeading tag="Stage 2 Output" tagClass="tag-amber"
                        title="Eᴵʲ Pairwise Distance Matrix" />
          <p className="text-xs mb-4 leading-5" style={{ color: "var(--text-secondary)" }}>
            Cα–Cα distances (Å) between pocket-lining residues in the selected transient conformation.
            Used as conditioning input for Stage 3 pharmacophore-guided SMILES generation.
          </p>
          <div className="flex gap-4 text-[11px] font-mono mb-4" style={{ color: "var(--text-muted)" }}>
            <span>Shape <span style={{ color: "var(--text-secondary)" }}>{(eij.shape ?? [0,0])[0]} × {(eij.shape ?? [0,0])[1]}</span></span>
            <span>Source <span style={{ color: "var(--text-secondary)" }}>{eij.file ?? "—"}</span></span>
          </div>
          <EijHeatmap labels={eij.residue_labels ?? []} values={eij.values ?? []} />
        </Panel>

        {/* ── Sequence ── */}
        <Panel>
          <PanelHeading tag="Sequence" tagClass="tag-green"
                        title="Full Protein Sequence" />
          <p className="font-mono text-[10px] mb-3" style={{ color: "var(--text-muted)" }}>
            {seq.fasta_header ?? `>${meta.gene ?? "?"} ${meta.mutation_mature ?? ""}`}
          </p>
          <div className="border p-4 overflow-x-auto"
               style={{ background: "var(--bg-1)", borderColor: "var(--border)", borderRadius: "2px" }}>
            <div className="font-mono text-[11px] break-all leading-5">
              {(seq.fasta_sequence ?? seq.full_sequence ?? "").split("").map((ch: string, i: number) => {
                const pos = i + 1;
                const isLining = seq.pocket_lining_positions_precursor?.includes(pos);
                return (
                  <span key={i}
                        style={{ color: isLining ? "var(--accent)" : "var(--text-muted)", fontWeight: isLining ? 700 : 400 }}
                        title={isLining ? `Pos ${pos} — pocket-lining` : undefined}>
                    {ch}
                  </span>
                );
              })}
            </div>
            <div className="mt-3 flex items-center gap-2">
              <span className="inline-block w-3 h-2" style={{ background: "var(--accent)" }} />
              <span className="text-[10px]" style={{ color: "var(--text-muted)" }}>
                Pocket-lining residues: <span className="font-mono font-semibold" style={{ color: "var(--text-secondary)" }}>
                  {seq.pocket_lining_1letter ?? "—"}
                </span>
              </span>
            </div>
          </div>
        </Panel>

        {/* ── Data Assets ── */}
        <Panel>
          <PanelHeading tag="Data Assets" tagClass="tag-cyan" title="Downloadable Files" />
          <div className="grid sm:grid-cols-2 gap-2">
            {[
              { label: "E_ij Matrix (.npy)", path: assets.eij_matrix_npy },
              ...(assets.eij_matrix_csv ? [{ label: "E_ij Matrix (.csv)", path: assets.eij_matrix_csv }] : []),
              { label: "Transient Conformation (.pdb)", path: assets.transient_conformation_pdb },
              { label: "Pocket Summary (.json)", path: assets.pocket_summary_json },
              { label: "Stage 3 Candidates (.json)", path: assets.stage3_candidates_json },
            ].map(a => (
              <div key={a.label} className="p-3 border" style={{ background: "var(--bg-1)", borderColor: "var(--border)", borderRadius: "2px" }}>
                <div className="data-label mb-1">{a.label}</div>
                <div className="font-mono text-[10px] break-all" style={{ color: "var(--text-secondary)" }}>{a.path}</div>
              </div>
            ))}
          </div>
        </Panel>

        {/* ── Footer nav ── */}
        <div className="flex justify-between items-center pt-1 pb-2">
          <Link href="/database"
                className="px-4 py-2 text-xs border font-medium transition-colors"
                style={{ borderColor: "var(--border-mid)", color: "var(--text-secondary)", borderRadius: "2px", background: "var(--bg)" }}>
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
