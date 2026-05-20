import Link from "next/link";

// ── Pipeline SVG diagram — fully self-contained, no external deps ────────────
function PipelineDiagram() {
  const box = (x: number, y: number, w: number, h: number, label: string, sub: string, color: string, tag?: string) => (
    <g key={label}>
      <rect x={x} y={y} width={w} height={h} rx={8}
            fill={color === "cyan" ? "#0c1f2e" : color === "violet" ? "#140e24" : "#0e1a0e"}
            stroke={color === "cyan" ? "#0ea5e9" : color === "violet" ? "#8b5cf6" : "#22c55e"}
            strokeWidth={1.5} />
      {tag && (
        <text x={x + 10} y={y + 14} fontSize={8} fontFamily="monospace"
              fill={color === "cyan" ? "#0ea5e9" : color === "violet" ? "#8b5cf6" : "#22c55e"}>
          {tag}
        </text>
      )}
      <text x={x + w / 2} y={y + (tag ? h / 2 + 6 : h / 2 - 4)} fontSize={12} fontWeight="bold"
            textAnchor="middle" fill="white" fontFamily="system-ui">
        {label}
      </text>
      <text x={x + w / 2} y={y + (tag ? h / 2 + 21 : h / 2 + 11)} fontSize={9}
            textAnchor="middle" fill="#94a3b8" fontFamily="system-ui">
        {sub}
      </text>
    </g>
  );

  const arrow = (x1: number, y1: number, x2: number, y2: number) => (
    <g key={`${x1}-${x2}`}>
      <defs>
        <marker id="arr" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto">
          <path d="M0,0 L0,6 L8,3 z" fill="#334155" />
        </marker>
      </defs>
      <line x1={x1} y1={y1} x2={x2} y2={y2} stroke="#334155" strokeWidth={1.5} markerEnd="url(#arr)" />
    </g>
  );

  const label = (x: number, y: number, txt: string) => (
    <text key={txt} x={x} y={y} fontSize={9} textAnchor="middle" fill="#475569" fontFamily="monospace">
      {txt}
    </text>
  );

  return (
    <svg viewBox="0 0 880 210" className="w-full" style={{ maxHeight: 210 }}>
      {/* Input */}
      <rect x={10} y={75} width={110} height={60} rx={6} fill="#0a0f1a" stroke="#334155" strokeWidth={1} strokeDasharray="4 2" />
      <text x={65} y={99} fontSize={10} textAnchor="middle" fill="#64748b" fontFamily="monospace">INPUT</text>
      <text x={65} y={113} fontSize={9} textAnchor="middle" fill="#94a3b8" fontFamily="system-ui">Pathogenic</text>
      <text x={65} y={124} fontSize={9} textAnchor="middle" fill="#94a3b8" fontFamily="system-ui">missense variant</text>

      {arrow(120, 105, 148, 105)}
      {label(134, 100, "")}

      {/* Stage 1 */}
      {box(150, 55, 170, 100, "Stage 1", "GNN-LM Classifier", "cyan", "RESCUE FILTER")}
      <text x={235} y={130} fontSize={8} textAnchor="middle" fill="#0ea5e9" fontFamily="monospace">Rescue score ≥ 0.70</text>

      {arrow(320, 105, 348, 105)}
      {label(334, 98, "amenable")}

      {/* Stage 2 */}
      {box(350, 45, 200, 120, "Stage 2", "ANM + fpocket", "violet", "POCKET DETECTION")}
      <text x={450} y={138} fontSize={8} textAnchor="middle" fill="#8b5cf6" fontFamily="monospace">20 conformations · E_ij matrix</text>
      <text x={450} y={148} fontSize={8} textAnchor="middle" fill="#8b5cf6" fontFamily="monospace">drug. threshold &gt; 0.70</text>

      {arrow(550, 105, 578, 105)}
      {label(564, 98, "cryptic pocket")}

      {/* Stage 3 */}
      {box(580, 55, 190, 100, "Stage 3", "SMILES Evolution", "green", "CHAPERONE DESIGN")}
      <text x={675} y={130} fontSize={8} textAnchor="middle" fill="#22c55e" fontFamily="monospace">MW &lt;350 Da · SA &lt;3.5 · QED &gt;0.7</text>

      {arrow(770, 105, 798, 105)}

      {/* Output */}
      <rect x={800} y={75} width={72} height={60} rx={6} fill="#0a1a0e" stroke="#22c55e" strokeWidth={1.5} />
      <text x={836} y={99} fontSize={9} textAnchor="middle" fill="#22c55e" fontFamily="monospace">OUTPUT</text>
      <text x={836} y={113} fontSize={8} textAnchor="middle" fill="#94a3b8" fontFamily="system-ui">Chaperone</text>
      <text x={836} y={124} fontSize={8} textAnchor="middle" fill="#94a3b8" fontFamily="system-ui">+ PDB atlas</text>

      {/* Stage labels at bottom */}
      <text x={235} y={200} fontSize={9} textAnchor="middle" fill="#475569" fontFamily="monospace">Multi-modal GNN ×</text>
      <text x={235} y={210} fontSize={9} textAnchor="middle" fill="#475569" fontFamily="monospace">Transformer fusion</text>
      <text x={450} y={200} fontSize={9} textAnchor="middle" fill="#475569" fontFamily="monospace">Anisotropic Network Model ·</text>
      <text x={450} y={210} fontSize={9} textAnchor="middle" fill="#475569" fontFamily="monospace">Voronoi alpha-sphere clustering</text>
      <text x={675} y={200} fontSize={9} textAnchor="middle" fill="#475569" fontFamily="monospace">Evolutionary SMILES search ·</text>
      <text x={675} y={210} fontSize={9} textAnchor="middle" fill="#475569" fontFamily="monospace">RDKit pharmacophore scoring</text>
    </svg>
  );
}

function StatPill({ label, value, color = "cyan" }: { label: string; value: string; color?: string }) {
  const c = color === "cyan" ? { border: "#0ea5e950", bg: "#0c1f2e", text: "#0ea5e9" }
           : color === "violet" ? { border: "#8b5cf650", bg: "#140e24", text: "#8b5cf6" }
           : { border: "#22c55e50", bg: "#0e1a0e", text: "#22c55e" };
  return (
    <div className="rounded-lg px-4 py-3 text-center border"
         style={{ background: c.bg, borderColor: c.border }}>
      <div className="text-xl font-black font-mono" style={{ color: c.text }}>{value}</div>
      <div className="text-[10px] mt-0.5" style={{ color: "#64748b", textTransform: "uppercase", letterSpacing: "0.08em" }}>{label}</div>
    </div>
  );
}

function ValidationCard({ gene, mutation, disease, drugScore, mw, sa, qed, pocket, residues, smiles }: {
  gene: string; mutation: string; disease: string; drugScore: string;
  mw: string; sa: string; qed: string; pocket: string; residues: string; smiles: string;
}) {
  return (
    <div className="rounded-xl border p-5 space-y-4" style={{ borderColor: "#1e293b", background: "#0a0f1a" }}>
      <div className="flex flex-wrap items-center gap-2">
        <span className="rounded px-2 py-0.5 text-xs font-bold font-mono border" style={{ background: "#0c1f2e", borderColor: "#0ea5e930", color: "#0ea5e9" }}>{gene}</span>
        <span className="rounded px-2 py-0.5 text-xs font-bold font-mono border" style={{ background: "#1c1208", borderColor: "#f59e0b30", color: "#f59e0b" }}>{mutation}</span>
        <span className="text-sm font-semibold" style={{ color: "#e2e8f0" }}>{disease}</span>
      </div>
      <div className="grid grid-cols-4 gap-2">
        {[
          { l: "Drug.", v: drugScore, c: "#0ea5e9" },
          { l: "MW (Da)", v: mw, c: "#94a3b8" },
          { l: "SA Score", v: sa, c: "#22c55e" },
          { l: "QED", v: qed, c: "#8b5cf6" },
        ].map(p => (
          <div key={p.l} className="rounded-lg p-2.5 text-center border" style={{ background: "#0d1117", borderColor: "#1e293b" }}>
            <div className="text-sm font-bold font-mono" style={{ color: p.c }}>{p.v}</div>
            <div className="text-[9px] mt-0.5" style={{ color: "#475569", textTransform: "uppercase" }}>{p.l}</div>
          </div>
        ))}
      </div>
      <div>
        <div className="text-[10px] mb-1" style={{ color: "#475569", textTransform: "uppercase", letterSpacing: "0.06em" }}>Cryptic Pocket</div>
        <div className="text-xs" style={{ color: "#64748b" }}>{pocket}</div>
      </div>
      <div>
        <div className="text-[10px] mb-1.5" style={{ color: "#475569", textTransform: "uppercase", letterSpacing: "0.06em" }}>Pocket-Lining Residues</div>
        <div className="flex flex-wrap gap-1">
          {residues.split("-").map(r => (
            <span key={r} className="rounded px-1.5 py-0.5 text-[10px] font-mono border" style={{ background: "#140e24", borderColor: "#8b5cf630", color: "#8b5cf6" }}>{r}</span>
          ))}
        </div>
      </div>
      <div>
        <div className="text-[10px] mb-1" style={{ color: "#475569", textTransform: "uppercase", letterSpacing: "0.06em" }}>Generated Chaperone (SMILES)</div>
        <div className="font-mono text-[10px] break-all leading-4 rounded p-2" style={{ background: "#0d1117", color: "#0ea5e9" }}>{smiles}</div>
      </div>
    </div>
  );
}

export default function AbstractPage() {
  return (
    <div style={{ background: "#07090f", minHeight: "100vh", color: "#e2e8f0" }}>

      {/* ── Top bar ── */}
      <div className="border-b" style={{ borderColor: "#1e293b" }}>
        <div className="max-w-4xl mx-auto px-6 h-12 flex items-center justify-between">
          <div className="flex items-center gap-3 text-xs font-mono" style={{ color: "#475569" }}>
            <span style={{ color: "#0ea5e9" }}>S.Y.A.L.I.S Labs</span>
            <span>·</span>
            <span>Preprint — 2026</span>
          </div>
          <div className="flex items-center gap-4 text-xs">
            <Link href="/database" style={{ color: "#475569" }} className="hover:text-white transition-colors">
              Database →
            </Link>
            <a href="https://github.com/aesenthilvanan-coder/REFOLD"
               target="_blank" rel="noopener noreferrer"
               className="flex items-center gap-1.5 px-3 py-1 rounded border text-xs transition-colors hover:border-white/20"
               style={{ borderColor: "#1e293b", color: "#94a3b8" }}>
              <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                <path d="M12 2C6.477 2 2 6.477 2 12c0 4.42 2.87 8.17 6.84 9.5.5.08.66-.23.66-.5v-1.69c-2.77.6-3.36-1.34-3.36-1.34-.46-1.16-1.11-1.47-1.11-1.47-.91-.62.07-.6.07-.6 1 .07 1.53 1.03 1.53 1.03.87 1.52 2.34 1.07 2.91.83.09-.65.35-1.09.63-1.34-2.22-.25-4.55-1.11-4.55-4.92 0-1.11.38-2 1.03-2.71-.1-.25-.45-1.29.1-2.64 0 0 .84-.27 2.75 1.02.79-.22 1.65-.33 2.5-.33.85 0 1.71.11 2.5.33 1.91-1.29 2.75-1.02 2.75-1.02.55 1.35.2 2.39.1 2.64.65.71 1.03 1.6 1.03 2.71 0 3.82-2.34 4.66-4.57 4.91.36.31.69.92.69 1.85V21c0 .27.16.59.67.5C19.14 20.16 22 16.42 22 12A10 10 0 0012 2z"/>
              </svg>
              GitHub
            </a>
          </div>
        </div>
      </div>

      <main className="max-w-4xl mx-auto px-6 py-14 space-y-14">

        {/* ── Title block ── */}
        <div className="space-y-6 text-center">
          <div className="space-y-1">
            <div className="text-xs font-mono tracking-widest uppercase" style={{ color: "#475569" }}>
              Computational Drug Discovery · Proteomics · ML
            </div>
          </div>

          <div>
            <div className="text-4xl font-black tracking-tight leading-tight" style={{ color: "#e2e8f0" }}>
              <span style={{ color: "#0ea5e9" }}>R.E.F.O.L.D</span>
            </div>
            <div className="text-base font-semibold mt-2 leading-snug" style={{ color: "#94a3b8" }}>
              Proteome-Scale De Novo Design of Pharmacological Chaperones<br />
              for Context-Locked Cryptic Pockets
            </div>
          </div>

          <div className="text-sm" style={{ color: "#475569" }}>
            Aaryan Senthilvanan &nbsp;·&nbsp; S.Y.A.L.I.S Labs &nbsp;·&nbsp; 2026
          </div>

          <div className="flex justify-center gap-3 flex-wrap">
            <a href="https://github.com/aesenthilvanan-coder/REFOLD"
               target="_blank" rel="noopener noreferrer"
               className="flex items-center gap-2 px-4 py-2 rounded-lg text-sm font-semibold border transition-all hover:border-sky-500/50"
               style={{ background: "#0c1f2e", borderColor: "#0ea5e930", color: "#0ea5e9" }}>
              <svg width="15" height="15" viewBox="0 0 24 24" fill="currentColor">
                <path d="M12 2C6.477 2 2 6.477 2 12c0 4.42 2.87 8.17 6.84 9.5.5.08.66-.23.66-.5v-1.69c-2.77.6-3.36-1.34-3.36-1.34-.46-1.16-1.11-1.47-1.11-1.47-.91-.62.07-.6.07-.6 1 .07 1.53 1.03 1.53 1.03.87 1.52 2.34 1.07 2.91.83.09-.65.35-1.09.63-1.34-2.22-.25-4.55-1.11-4.55-4.92 0-1.11.38-2 1.03-2.71-.1-.25-.45-1.29.1-2.64 0 0 .84-.27 2.75 1.02.79-.22 1.65-.33 2.5-.33.85 0 1.71.11 2.5.33 1.91-1.29 2.75-1.02 2.75-1.02.55 1.35.2 2.39.1 2.64.65.71 1.03 1.6 1.03 2.71 0 3.82-2.34 4.66-4.57 4.91.36.31.69.92.69 1.85V21c0 .27.16.59.67.5C19.14 20.16 22 16.42 22 12A10 10 0 0012 2z"/>
              </svg>
              View Source — github.com/aesenthilvanan-coder/REFOLD
            </a>
            <Link href="/database"
               className="flex items-center gap-2 px-4 py-2 rounded-lg text-sm font-semibold border transition-all hover:border-violet-500/50"
               style={{ background: "#140e24", borderColor: "#8b5cf630", color: "#8b5cf6" }}>
              Browse PCD Database →
            </Link>
          </div>
        </div>

        {/* ── Live stats strip ── */}
        <div className="grid grid-cols-2 sm:grid-cols-4 gap-3">
          <StatPill label="Database Entries" value="668+" color="cyan" />
          <StatPill label="Diseases Covered" value="231+" color="violet" />
          <StatPill label="Avg Drug. Score" value="0.835" color="cyan" />
          <StatPill label="ClinVar Variants" value="178,597" color="violet" />
        </div>

        {/* ── Abstract ── */}
        <div className="space-y-4">
          <h2 className="text-xs font-mono tracking-widest uppercase" style={{ color: "#475569" }}>Abstract</h2>
          <div className="space-y-4 text-sm leading-7" style={{ color: "#94a3b8" }}>
            <p>
              Pathogenic missense mutations frequently cause disease through accelerated protein misfolding rather than direct functional disruption — a mechanism that leaves the vast majority of rare disease variants without any approved pharmacological intervention. Traditional structure-based drug discovery remains largely blind to these variants, relying on static wild-type templates that fail to capture the transient conformational trajectories unique to destabilized mutant ensembles.
            </p>
            <p>
              To overcome this, we introduce <span style={{ color: "#e2e8f0", fontWeight: 600 }}>REFOLD</span>, a fully automated computational framework for the proteome-scale <em>de novo</em> synthesis of pharmacological chaperones. REFOLD couples a <span style={{ color: "#0ea5e9" }}>multi-modal GNN-Transformer fusion classifier</span> for rescue amenability scoring with an <span style={{ color: "#8b5cf6" }}>Anisotropic Network Model (ANM)</span> that simulates 20 mutant transition-state conformations at 1.5 Å RMSD to expose hidden, allosteric cryptic pockets undetectable in wild-type structures. A <span style={{ color: "#22c55e" }}>pharmacophore-guided evolutionary SMILES search</span>, scored against RDKit molecular descriptors and the E<sub>ij</sub> pairwise Cα–Cα distance matrix of pocket-lining residues, then designs highly druglike, synthetically accessible small molecules tailored to these coordinate constraints entirely from scratch.
            </p>
            <p>
              We validate REFOLD on two historically intractable Class II misfolding targets. For the <span style={{ color: "#e2e8f0" }}>GBA1 L444P variant</span> (Gaucher disease Type I), the model bypassed the collapsed canonical active site (WT fpocket druggability: 0.11), identified a distal allosteric hinge at residues A40–R41–P42–C43–D63–S64–F65–R87–M88–E89–L90 (mutant pocket druggability: <span style={{ color: "#0ea5e9" }}>0.926</span>, volume: 462.6 Å³), and generated a bis-aromatic piperidine scaffold with MW 317 Da, SA score 2.7, and QED 0.79. For the <span style={{ color: "#e2e8f0" }}>CFTR G85E variant</span> (Cystic Fibrosis), REFOLD similarly bypassed the VX-809 canonical binding site (druggability 0.003 in mutant), instead targeting the exposed ICL1-TM2 junction (druggability: <span style={{ color: "#0ea5e9" }}>0.807</span>, volume: 484.0 Å³) with a fluorinated piperazine scaffold of MW 316 Da and exceptional SA score <span style={{ color: "#22c55e" }}>1.7</span>. In both cases, molecules act as physical &ldquo;kinetic splints&rdquo; that stabilize misfolded intermediates and rescue structural integrity before ER quality control-mediated degradation.
            </p>
            <p>
              Scaled across the full ClinVar pathogenic missense catalog (<span style={{ color: "#0ea5e9" }}>178,597 variants</span>), REFOLD has to date generated <span style={{ color: "#0ea5e9" }}>668 complete chaperone entries</span> spanning <span style={{ color: "#8b5cf6" }}>231 distinct diseases</span>, with a mean pocket druggability of 0.835 across all accepted variants. All results are continuously published to the open-access <span style={{ color: "#e2e8f0" }}>Pharmacological Chaperone Database (PCD)</span>, providing interactive 3D transient-conformation visualizations, full E<sub>ij</sub> matrices, molecular property profiles, and SMILES strings for every entry — establishing a scalable blueprint for zero-shot orphan disease drug discovery.
            </p>
          </div>
        </div>

        {/* ── Architecture ── */}
        <div className="space-y-5">
          <h2 className="text-xs font-mono tracking-widest uppercase" style={{ color: "#475569" }}>Pipeline Architecture</h2>

          <div className="rounded-xl border p-5" style={{ borderColor: "#1e293b", background: "#0a0f1a" }}>
            <PipelineDiagram />
          </div>

          <div className="grid sm:grid-cols-3 gap-4 text-sm">
            {[
              {
                stage: "Stage 1",
                tag: "Rescue Filter",
                color: "#0ea5e9",
                bg: "#0c1f2e",
                border: "#0ea5e930",
                title: "GNN-LM Fusion Classifier",
                body: "A multi-modal graph neural network encodes the local residue contact graph of the mutant structure while a Transformer-based language model captures sequence-level evolutionary context. The fused representation produces a rescue amenability score; only variants with score ≥ 0.70 proceed. This eliminates loss-of-function variants that chaperones cannot rescue.",
              },
              {
                stage: "Stage 2",
                tag: "Pocket Detection",
                color: "#8b5cf6",
                bg: "#140e24",
                border: "#8b5cf630",
                title: "ANM + fpocket",
                body: "An Anisotropic Network Model samples 20 mutant transition-state conformations along low-frequency eigenmodes at 1.5 Å RMSD from the AlphaFold structure. fpocket runs Voronoi alpha-sphere clustering on each conformation; pockets absent from the WT ensemble but present in ≥ 3 mutant conformations are flagged as cryptic. The E_ij pairwise Cα–Cα distance matrix of lining residues is extracted as a pharmacophore conditioning signal.",
              },
              {
                stage: "Stage 3",
                tag: "Chaperone Design",
                color: "#22c55e",
                bg: "#0e1a0e",
                border: "#22c55e30",
                title: "Evolutionary SMILES Search",
                body: "Starting from a diverse seed population, an evolutionary SMILES search applies mutation, crossover, and elite selection over multiple generations. Each candidate is scored on a composite of fpocket affinity estimate, Lipinski Ro5, Veber oral bioavailability, SA score (< 3.5), QED (> 0.7), and pharmacophore complementarity to the E_ij matrix. Top-10 candidates per variant are stored.",
              },
            ].map(s => (
              <div key={s.stage} className="rounded-xl border p-4 space-y-2"
                   style={{ background: s.bg, borderColor: s.border }}>
                <div className="flex items-center gap-2">
                  <span className="text-[9px] font-mono px-2 py-0.5 rounded border"
                        style={{ background: "#0d1117", borderColor: s.border, color: s.color }}>
                    {s.stage} · {s.tag}
                  </span>
                </div>
                <div className="text-sm font-semibold" style={{ color: "#e2e8f0" }}>{s.title}</div>
                <div className="text-xs leading-5" style={{ color: "#64748b" }}>{s.body}</div>
              </div>
            ))}
          </div>
        </div>

        {/* ── Validation ── */}
        <div className="space-y-5">
          <h2 className="text-xs font-mono tracking-widest uppercase" style={{ color: "#475569" }}>Validation Cases</h2>
          <div className="grid sm:grid-cols-2 gap-5">
            <ValidationCard
              gene="GBA1" mutation="L444P"
              disease="Gaucher Disease Type I"
              drugScore="0.926" mw="317" sa="2.7" qed="0.795"
              pocket="Distal allosteric hinge — absent in WT (drug=0.11), exposed by L444P ensemble"
              residues="A40-R41-P42-C43-D63-S64-F65-R87-M88-E89-L90"
              smiles="O=C(NCc1ccc(O)c(N2CCCCC2)c1)C1CCNCC1"
            />
            <ValidationCard
              gene="CFTR" mutation="G85E"
              disease="Cystic Fibrosis"
              drugScore="0.807" mw="316" sa="1.7" qed="0.869"
              pocket="ICL1-TM2 junction — VX-809 canonical site destroyed (drug=0.003)"
              residues="I70-L73-R74-F77-F78-F81-L195"
              smiles="O=C(c1cc(F)ccc1F)N1CCN(Cc2ccccc2)CC1"
            />
          </div>
        </div>

        {/* ── Key design principles ── */}
        <div className="space-y-4">
          <h2 className="text-xs font-mono tracking-widest uppercase" style={{ color: "#475569" }}>Key Design Principles</h2>
          <div className="grid sm:grid-cols-2 gap-3 text-xs leading-6" style={{ color: "#64748b" }}>
            {[
              ["Context-locked cryptic pockets", "Pockets are only sampled in mutant conformational trajectories — not visible in WT crystal structures or static AlphaFold models. REFOLD specifically targets this structurally dynamic regime."],
              ["Kinetic splint mechanism", "Generated chaperones are not designed to restore enzymatic function directly; they stabilize misfolded intermediates long enough to evade ERAD and traffic correctly to their destination organelle."],
              ["Fully de novo, zero-shot", "No experimental binding data, crystal co-complex, or known ligand scaffold is used. Every chaperone is generated from scratch conditioned solely on the cryptic pocket geometry."],
              ["Proteome-scale automation", "The pipeline runs unattended: ClinVar variant queue → AlphaFold fetch → ANM sampling → fpocket screening → SMILES evolution → GitHub push → live website update. One full entry takes ~2 minutes."],
            ].map(([title, body]) => (
              <div key={title as string} className="rounded-xl border p-4" style={{ borderColor: "#1e293b", background: "#0a0f1a" }}>
                <div className="text-xs font-semibold mb-1" style={{ color: "#e2e8f0" }}>{title as string}</div>
                <div>{body as string}</div>
              </div>
            ))}
          </div>
        </div>

        {/* ── Footer ── */}
        <div className="border-t pt-8 flex flex-col sm:flex-row justify-between items-start gap-4"
             style={{ borderColor: "#1e293b" }}>
          <div className="text-xs space-y-1" style={{ color: "#475569" }}>
            <div style={{ color: "#94a3b8" }}>Aaryan Senthilvanan · S.Y.A.L.I.S Labs · 2026</div>
            <div>All data, code, and generated structures are openly available.</div>
          </div>
          <div className="flex gap-3">
            <a href="https://github.com/aesenthilvanan-coder/REFOLD"
               target="_blank" rel="noopener noreferrer"
               className="px-4 py-2 rounded-lg text-xs font-semibold border transition-colors hover:border-sky-500/40"
               style={{ background: "#0c1f2e", borderColor: "#0ea5e930", color: "#0ea5e9" }}>
              GitHub →
            </a>
            <Link href="/database"
               className="px-4 py-2 rounded-lg text-xs font-semibold border transition-colors"
               style={{ background: "#0d1117", borderColor: "#1e293b", color: "#64748b" }}>
              PCD Database →
            </Link>
          </div>
        </div>

      </main>
    </div>
  );
}
