import Link from "next/link";

// ── Pipeline SVG — light-mode institutional palette ──────────────────────────
function PipelineDiagram() {
  const CY = 120;
  return (
    <svg viewBox="0 0 900 255" className="w-full" style={{ maxHeight: 240 }}>
      <defs>
        <marker id="arr" markerWidth="8" markerHeight="8" refX="6" refY="3" orient="auto">
          <path d="M0,0 L0,6 L8,3 z" fill="#94A3B8" />
        </marker>
      </defs>

      {/* Input */}
      <rect x={8} y={88} width={108} height={64} rx={3}
            fill="#F8FAFC" stroke="#CBD5E1" strokeWidth={1} strokeDasharray="4 2" />
      <text x={62} y={110} fontSize={9} textAnchor="middle" fill="#94A3B8" fontFamily="monospace" fontWeight="bold">INPUT</text>
      <text x={62} y={125} fontSize={9} textAnchor="middle" fill="#475569" fontFamily="system-ui">Pathogenic</text>
      <text x={62} y={138} fontSize={9} textAnchor="middle" fill="#475569" fontFamily="system-ui">missense variant</text>

      <line x1={116} y1={CY} x2={148} y2={CY} stroke="#94A3B8" strokeWidth={1.5} markerEnd="url(#arr)" />

      {/* Stage 1 */}
      <rect x={150} y={42} width={162} height={138} rx={3}
            fill="#EFF6FF" stroke="#1E3A8A" strokeWidth={1.5} />
      <text x={162} y={58} fontSize={8} fontFamily="monospace" fill="#1E3A8A" fontWeight="bold">RESCUE FILTER</text>
      <text x={231} y={98} fontSize={12} fontWeight="bold" textAnchor="middle" fill="#0F172A" fontFamily="system-ui">Stage 1</text>
      <text x={231} y={114} fontSize={10} textAnchor="middle" fill="#334155" fontFamily="system-ui">GNN-LM Classifier</text>
      <text x={231} y={130} fontSize={8} textAnchor="middle" fill="#1E3A8A" fontFamily="monospace">Rescue score ≥ 0.70</text>
      <text x={231} y={146} fontSize={8} textAnchor="middle" fill="#64748B" fontFamily="monospace">Multi-modal GNN × Transformer</text>
      <text x={231} y={158} fontSize={8} textAnchor="middle" fill="#64748B" fontFamily="monospace">sequence + structure fusion</text>

      <line x1={312} y1={CY} x2={346} y2={CY} stroke="#94A3B8" strokeWidth={1.5} markerEnd="url(#arr)" />
      <text x={329} y={113} fontSize={8} textAnchor="middle" fill="#94A3B8" fontFamily="monospace">amenable</text>

      {/* Stage 2 */}
      <rect x={348} y={28} width={200} height={166} rx={3}
            fill="#F0FDFA" stroke="#0F766E" strokeWidth={1.5} />
      <text x={360} y={44} fontSize={8} fontFamily="monospace" fill="#0F766E" fontWeight="bold">POCKET DETECTION</text>
      <text x={448} y={88} fontSize={12} fontWeight="bold" textAnchor="middle" fill="#0F172A" fontFamily="system-ui">Stage 2</text>
      <text x={448} y={104} fontSize={10} textAnchor="middle" fill="#334155" fontFamily="system-ui">ANM + fpocket</text>
      <text x={448} y={120} fontSize={8} textAnchor="middle" fill="#0F766E" fontFamily="monospace">20 conformations · 1.5 Å RMSD</text>
      <text x={448} y={134} fontSize={8} textAnchor="middle" fill="#0F766E" fontFamily="monospace">drug. threshold &gt; 0.70</text>
      <text x={448} y={151} fontSize={8} textAnchor="middle" fill="#64748B" fontFamily="monospace">Anisotropic Network Model ·</text>
      <text x={448} y={163} fontSize={8} textAnchor="middle" fill="#64748B" fontFamily="monospace">Voronoi alpha-sphere clustering</text>

      <line x1={548} y1={CY} x2={582} y2={CY} stroke="#94A3B8" strokeWidth={1.5} markerEnd="url(#arr)" />
      <text x={565} y={113} fontSize={8} textAnchor="middle" fill="#94A3B8" fontFamily="monospace">cryptic pocket</text>

      {/* Stage 3 */}
      <rect x={584} y={42} width={178} height={138} rx={3}
            fill="#F0FDF4" stroke="#15803D" strokeWidth={1.5} />
      <text x={596} y={58} fontSize={8} fontFamily="monospace" fill="#15803D" fontWeight="bold">CHAPERONE DESIGN</text>
      <text x={673} y={98} fontSize={12} fontWeight="bold" textAnchor="middle" fill="#0F172A" fontFamily="system-ui">Stage 3</text>
      <text x={673} y={114} fontSize={10} textAnchor="middle" fill="#334155" fontFamily="system-ui">SMILES Evolution</text>
      <text x={673} y={130} fontSize={8} textAnchor="middle" fill="#15803D" fontFamily="monospace">MW &lt;350 Da · SA &lt;3.5 · QED &gt;0.7</text>
      <text x={673} y={146} fontSize={8} textAnchor="middle" fill="#64748B" fontFamily="monospace">Evolutionary search ·</text>
      <text x={673} y={158} fontSize={8} textAnchor="middle" fill="#64748B" fontFamily="monospace">RDKit pharmacophore scoring</text>

      <line x1={762} y1={CY} x2={796} y2={CY} stroke="#94A3B8" strokeWidth={1.5} markerEnd="url(#arr)" />

      {/* Output */}
      <rect x={798} y={88} width={94} height={64} rx={3}
            fill="#F0FDF4" stroke="#15803D" strokeWidth={1.5} />
      <text x={845} y={110} fontSize={9} textAnchor="middle" fill="#15803D" fontFamily="monospace" fontWeight="bold">OUTPUT</text>
      <text x={845} y={125} fontSize={9} textAnchor="middle" fill="#334155" fontFamily="system-ui">Chaperone</text>
      <text x={845} y={138} fontSize={9} textAnchor="middle" fill="#334155" fontFamily="system-ui">+ PDB atlas</text>
    </svg>
  );
}

function StatPill({ label, value, color = "navy" }: { label: string; value: string; color?: "navy" | "teal" | "green" }) {
  const styles = {
    navy:  { border: "#BFDBFE", bg: "#EFF6FF", text: "#1E3A8A" },
    teal:  { border: "#99F6E4", bg: "#F0FDFA", text: "#0F766E" },
    green: { border: "#BBF7D0", bg: "#F0FDF4", text: "#15803D" },
  }[color];
  return (
    <div className="px-4 py-3 text-center border"
         style={{ background: styles.bg, borderColor: styles.border, borderRadius: "3px" }}>
      <div className="text-xl font-bold font-mono tabular-nums" style={{ color: styles.text }}>{value}</div>
      <div className="text-[10px] mt-0.5 font-semibold uppercase tracking-widest" style={{ color: "#64748B" }}>{label}</div>
    </div>
  );
}

function ValidationCard({ gene, mutation, disease, drugScore, mw, sa, qed, pocket, residues, smiles }: {
  gene: string; mutation: string; disease: string; drugScore: string;
  mw: string; sa: string; qed: string; pocket: string; residues: string; smiles: string;
}) {
  return (
    <div className="border p-5 space-y-4"
         style={{ borderColor: "#E2E8F0", background: "#FFFFFF", borderRadius: "3px" }}>
      <div className="flex flex-wrap items-center gap-2">
        <span className="tag tag-cyan">{gene}</span>
        <span className="tag tag-amber">{mutation}</span>
        <span className="text-sm font-semibold" style={{ color: "#0F172A" }}>{disease}</span>
      </div>

      <div className="grid grid-cols-4 gap-1.5">
        {[
          { l: "Drug.", v: drugScore, c: "#1E3A8A" },
          { l: "MW (Da)", v: mw,      c: "#334155" },
          { l: "SA Score", v: sa,     c: "#15803D" },
          { l: "QED",      v: qed,    c: "#0F766E" },
        ].map(p => (
          <div key={p.l} className="p-2.5 text-center border"
               style={{ background: "#F8FAFC", borderColor: "#E2E8F0", borderRadius: "2px" }}>
            <div className="text-sm font-bold font-mono tabular-nums" style={{ color: p.c }}>{p.v}</div>
            <div className="text-[9px] mt-0.5 font-semibold uppercase tracking-wider" style={{ color: "#94A3B8" }}>{p.l}</div>
          </div>
        ))}
      </div>

      <div>
        <div className="data-label mb-1">Cryptic Pocket</div>
        <div className="text-xs" style={{ color: "#475569" }}>{pocket}</div>
      </div>

      <div>
        <div className="data-label mb-2">Pocket-Lining Residues</div>
        <div className="flex flex-wrap gap-1">
          {residues.split("-").map(r => (
            <span key={r} className="tag tag-violet">{r}</span>
          ))}
        </div>
      </div>

      <div>
        <div className="data-label mb-1">Generated Chaperone (SMILES)</div>
        <div className="font-mono text-[10px] break-all leading-4 p-2 border"
             style={{ background: "#F8FAFC", borderColor: "#E2E8F0", color: "#1E3A8A", borderRadius: "2px" }}>
          {smiles}
        </div>
      </div>
    </div>
  );
}

export default function AbstractPage() {
  return (
    <div style={{ background: "#FFFFFF", minHeight: "100vh", color: "#0F172A" }}>

      {/* ── Top bar ── */}
      <div className="border-b" style={{ borderColor: "#E2E8F0" }}>
        <div className="max-w-4xl mx-auto px-6 h-12 flex items-center justify-between">
          <div className="flex items-center gap-3 text-xs font-mono" style={{ color: "#64748B" }}>
            <span className="font-bold" style={{ color: "#1E3A8A" }}>S.Y.A.L.I.S Labs</span>
            <span>·</span>
            <span>Preprint — 2026</span>
          </div>
          <div className="flex items-center gap-4 text-xs">
            <Link href="/database" style={{ color: "#64748B" }} className="hover:text-slate-900 transition-colors">
              Database →
            </Link>
            <a href="https://github.com/aesenthilvanan-coder/REFOLD"
               target="_blank" rel="noopener noreferrer"
               className="flex items-center gap-1.5 px-3 py-1 border text-xs transition-colors"
               style={{ borderColor: "#E2E8F0", color: "#475569", borderRadius: "2px" }}>
              <svg width="13" height="13" viewBox="0 0 24 24" fill="currentColor">
                <path d="M12 2C6.477 2 2 6.477 2 12c0 4.42 2.87 8.17 6.84 9.5.5.08.66-.23.66-.5v-1.69c-2.77.6-3.36-1.34-3.36-1.34-.46-1.16-1.11-1.47-1.11-1.47-.91-.62.07-.6.07-.6 1 .07 1.53 1.03 1.53 1.03.87 1.52 2.34 1.07 2.91.83.09-.65.35-1.09.63-1.34-2.22-.25-4.55-1.11-4.55-4.92 0-1.11.38-2 1.03-2.71-.1-.25-.45-1.29.1-2.64 0 0 .84-.27 2.75 1.02.79-.22 1.65-.33 2.5-.33.85 0 1.71.11 2.5.33 1.91-1.29 2.75-1.02 2.75-1.02.55 1.35.2 2.39.1 2.64.65.71 1.03 1.6 1.03 2.71 0 3.82-2.34 4.66-4.57 4.91.36.31.69.92.69 1.85V21c0 .27.16.59.67.5C19.14 20.16 22 16.42 22 12A10 10 0 0012 2z"/>
              </svg>
              GitHub
            </a>
          </div>
        </div>
      </div>

      <main className="max-w-4xl mx-auto px-6 py-12 space-y-12">

        {/* ── Title block ── */}
        <div className="space-y-5 text-center border-b pb-10" style={{ borderColor: "#E2E8F0" }}>
          <div className="text-[10px] font-mono tracking-widest uppercase" style={{ color: "#94A3B8" }}>
            Computational Drug Discovery · Proteomics · Machine Learning
          </div>
          <div>
            <div className="text-3xl font-bold tracking-tight" style={{ color: "#1E3A8A" }}>
              R.E.F.O.L.D
            </div>
            <div className="text-base font-semibold mt-2 leading-snug" style={{ color: "#334155" }}>
              Proteome-Scale De Novo Design of Pharmacological Chaperones<br />
              for Context-Locked Cryptic Pockets
            </div>
          </div>
          <div className="text-sm" style={{ color: "#64748B" }}>
            Aaryan Senthilvanan &nbsp;·&nbsp; S.Y.A.L.I.S Labs &nbsp;·&nbsp; 2026
          </div>
          <div className="flex justify-center gap-2 flex-wrap">
            <a href="https://github.com/aesenthilvanan-coder/REFOLD"
               target="_blank" rel="noopener noreferrer"
               className="flex items-center gap-2 px-4 py-2 border text-xs font-semibold transition-colors"
               style={{ background: "#EFF6FF", borderColor: "#BFDBFE", color: "#1E3A8A", borderRadius: "2px" }}>
              <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                <path d="M12 2C6.477 2 2 6.477 2 12c0 4.42 2.87 8.17 6.84 9.5.5.08.66-.23.66-.5v-1.69c-2.77.6-3.36-1.34-3.36-1.34-.46-1.16-1.11-1.47-1.11-1.47-.91-.62.07-.6.07-.6 1 .07 1.53 1.03 1.53 1.03.87 1.52 2.34 1.07 2.91.83.09-.65.35-1.09.63-1.34-2.22-.25-4.55-1.11-4.55-4.92 0-1.11.38-2 1.03-2.71-.1-.25-.45-1.29.1-2.64 0 0 .84-.27 2.75 1.02.79-.22 1.65-.33 2.5-.33.85 0 1.71.11 2.5.33 1.91-1.29 2.75-1.02 2.75-1.02.55 1.35.2 2.39.1 2.64.65.71 1.03 1.6 1.03 2.71 0 3.82-2.34 4.66-4.57 4.91.36.31.69.92.69 1.85V21c0 .27.16.59.67.5C19.14 20.16 22 16.42 22 12A10 10 0 0012 2z"/>
              </svg>
              View Source — github.com/aesenthilvanan-coder/REFOLD
            </a>
            <Link href="/database"
               className="flex items-center gap-2 px-4 py-2 border text-xs font-semibold transition-colors"
               style={{ background: "#F0FDFA", borderColor: "#99F6E4", color: "#0F766E", borderRadius: "2px" }}>
              Browse PCD Database →
            </Link>
          </div>
        </div>

        {/* ── Live stats strip ── */}
        <div className="grid grid-cols-2 sm:grid-cols-4 gap-2">
          <StatPill label="Database Entries" value="668+" color="navy" />
          <StatPill label="Diseases Covered" value="231+" color="teal" />
          <StatPill label="Avg Drug. Score" value="0.835" color="navy" />
          <StatPill label="ClinVar Variants" value="178,597" color="teal" />
        </div>

        {/* ── Abstract ── */}
        <section className="space-y-4">
          <h2 className="text-[10px] font-mono tracking-widest uppercase border-b pb-2"
              style={{ color: "#94A3B8", borderColor: "#E2E8F0" }}>Abstract</h2>
          <div className="space-y-3 text-sm leading-7" style={{ color: "#334155" }}>
            <p>
              Pathogenic missense mutations frequently cause disease through accelerated protein misfolding rather than direct functional disruption — a mechanism that leaves the vast majority of rare disease variants without any approved pharmacological intervention. Traditional structure-based drug discovery remains largely blind to these variants, relying on static wild-type templates that fail to capture the transient conformational trajectories unique to destabilized mutant ensembles.
            </p>
            <p>
              To overcome this, we introduce <span style={{ color: "#0F172A", fontWeight: 600 }}>REFOLD</span>, a fully automated computational framework for the proteome-scale <em>de novo</em> synthesis of pharmacological chaperones. REFOLD couples a <span style={{ color: "#1E3A8A", fontWeight: 600 }}>multi-modal GNN-Transformer fusion classifier</span> for rescue amenability scoring with an <span style={{ color: "#0F766E", fontWeight: 600 }}>Anisotropic Network Model (ANM)</span> that simulates 20 mutant transition-state conformations at 1.5 Å RMSD to expose hidden, allosteric cryptic pockets undetectable in wild-type structures. A <span style={{ color: "#15803D", fontWeight: 600 }}>pharmacophore-guided evolutionary SMILES search</span>, scored against RDKit molecular descriptors and the E<sub>ij</sub> pairwise Cα–Cα distance matrix of pocket-lining residues, then designs highly druglike, synthetically accessible small molecules tailored to these coordinate constraints entirely from scratch.
            </p>
            <p>
              We validate REFOLD on two historically intractable Class II misfolding targets. For the <span style={{ color: "#0F172A", fontWeight: 600 }}>GBA1 L444P variant</span> (Gaucher disease Type I), the model bypassed the collapsed canonical active site (WT fpocket druggability: 0.11), identified a distal allosteric hinge at residues A40–R41–P42–C43–D63–S64–F65–R87–M88–E89–L90 (mutant pocket druggability: <span style={{ color: "#1E3A8A", fontWeight: 600 }}>0.926</span>, volume: 462.6 Å³), and generated a bis-aromatic piperidine scaffold with MW 317 Da, SA score 2.7, and QED 0.79. For the <span style={{ color: "#0F172A", fontWeight: 600 }}>CFTR G85E variant</span> (Cystic Fibrosis), REFOLD similarly bypassed the VX-809 canonical binding site (druggability 0.003 in mutant), instead targeting the exposed ICL1-TM2 junction (druggability: <span style={{ color: "#1E3A8A", fontWeight: 600 }}>0.807</span>, volume: 484.0 Å³) with a fluorinated piperazine scaffold of MW 316 Da and exceptional SA score <span style={{ color: "#15803D", fontWeight: 600 }}>1.7</span>. In both cases, molecules act as physical &ldquo;kinetic splints&rdquo; that stabilize misfolded intermediates and rescue structural integrity before ER quality control-mediated degradation.
            </p>
            <p>
              Scaled across the full ClinVar pathogenic missense catalog (<span style={{ color: "#1E3A8A", fontWeight: 600 }}>178,597 variants</span>), REFOLD has to date generated <span style={{ color: "#1E3A8A", fontWeight: 600 }}>668 complete chaperone entries</span> spanning <span style={{ color: "#0F766E", fontWeight: 600 }}>231 distinct diseases</span>, with a mean pocket druggability of 0.835 across all accepted variants. All results are continuously published to the open-access <span style={{ color: "#0F172A", fontWeight: 600 }}>Pharmacological Chaperone Database (PCD)</span>, providing interactive 3D transient-conformation visualizations, full E<sub>ij</sub> matrices, molecular property profiles, and SMILES strings for every entry — establishing a scalable blueprint for zero-shot orphan disease drug discovery.
            </p>
          </div>
        </section>

        {/* ── Architecture ── */}
        <section className="space-y-5">
          <h2 className="text-[10px] font-mono tracking-widest uppercase border-b pb-2"
              style={{ color: "#94A3B8", borderColor: "#E2E8F0" }}>Pipeline Architecture</h2>

          <div className="border p-5" style={{ borderColor: "#E2E8F0", background: "#F8FAFC", borderRadius: "3px" }}>
            <PipelineDiagram />
          </div>

          <div className="grid sm:grid-cols-3 gap-3 text-sm">
            {[
              {
                stage: "Stage 1", tag: "Rescue Filter",
                border: "#BFDBFE", bg: "#EFF6FF", tagColor: "#1E3A8A",
                title: "GNN-LM Fusion Classifier",
                body: "A multi-modal graph neural network encodes the local residue contact graph of the mutant structure while a Transformer-based language model captures sequence-level evolutionary context. The fused representation produces a rescue amenability score; only variants with score ≥ 0.70 proceed. This eliminates loss-of-function variants that chaperones cannot rescue.",
              },
              {
                stage: "Stage 2", tag: "Pocket Detection",
                border: "#99F6E4", bg: "#F0FDFA", tagColor: "#0F766E",
                title: "ANM + fpocket",
                body: "An Anisotropic Network Model samples 20 mutant transition-state conformations along low-frequency eigenmodes at 1.5 Å RMSD from the AlphaFold structure. fpocket runs Voronoi alpha-sphere clustering on each conformation; pockets absent from the WT ensemble but present in ≥ 3 mutant conformations are flagged as cryptic. The E_ij pairwise Cα–Cα distance matrix of lining residues is extracted as a pharmacophore conditioning signal.",
              },
              {
                stage: "Stage 3", tag: "Chaperone Design",
                border: "#BBF7D0", bg: "#F0FDF4", tagColor: "#15803D",
                title: "Evolutionary SMILES Search",
                body: "Starting from a diverse seed population, an evolutionary SMILES search applies mutation, crossover, and elite selection over multiple generations. Each candidate is scored on a composite of fpocket affinity estimate, Lipinski Ro5, Veber oral bioavailability, SA score (< 3.5), QED (> 0.7), and pharmacophore complementarity to the E_ij matrix. Top-10 candidates per variant are stored.",
              },
            ].map(s => (
              <div key={s.stage} className="border p-4 space-y-2"
                   style={{ background: s.bg, borderColor: s.border, borderRadius: "3px" }}>
                <div className="flex items-center gap-2">
                  <span className="tag text-[9px]"
                        style={{ background: "#FFFFFF", borderColor: s.border, color: s.tagColor }}>
                    {s.stage} · {s.tag}
                  </span>
                </div>
                <div className="text-sm font-semibold" style={{ color: "#0F172A" }}>{s.title}</div>
                <div className="text-xs leading-5" style={{ color: "#475569" }}>{s.body}</div>
              </div>
            ))}
          </div>
        </section>

        {/* ── Validation cases ── */}
        <section className="space-y-5">
          <h2 className="text-[10px] font-mono tracking-widest uppercase border-b pb-2"
              style={{ color: "#94A3B8", borderColor: "#E2E8F0" }}>Validation Cases</h2>
          <div className="grid sm:grid-cols-2 gap-4">
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
        </section>

        {/* ── Key design principles ── */}
        <section className="space-y-4">
          <h2 className="text-[10px] font-mono tracking-widest uppercase border-b pb-2"
              style={{ color: "#94A3B8", borderColor: "#E2E8F0" }}>Key Design Principles</h2>
          <div className="grid sm:grid-cols-2 gap-3">
            {[
              ["Context-locked cryptic pockets", "Pockets are only sampled in mutant conformational trajectories — not visible in WT crystal structures or static AlphaFold models. REFOLD specifically targets this structurally dynamic regime."],
              ["Kinetic splint mechanism", "Generated chaperones are not designed to restore enzymatic function directly; they stabilize misfolded intermediates long enough to evade ERAD and traffic correctly to their destination organelle."],
              ["Fully de novo, zero-shot", "No experimental binding data, crystal co-complex, or known ligand scaffold is used. Every chaperone is generated from scratch conditioned solely on the cryptic pocket geometry."],
              ["Proteome-scale automation", "The pipeline runs unattended: ClinVar variant queue → AlphaFold fetch → ANM sampling → fpocket screening → SMILES evolution → GitHub push → live website update. One full entry takes ~2 minutes."],
            ].map(([title, body]) => (
              <div key={title as string} className="border p-4"
                   style={{ borderColor: "#E2E8F0", background: "#F8FAFC", borderRadius: "3px" }}>
                <div className="text-xs font-semibold mb-1.5" style={{ color: "#0F172A" }}>{title as string}</div>
                <div className="text-xs leading-5" style={{ color: "#475569" }}>{body as string}</div>
              </div>
            ))}
          </div>
        </section>

        {/* ── Footer ── */}
        <div className="border-t pt-6 flex flex-col sm:flex-row justify-between items-start gap-4"
             style={{ borderColor: "#E2E8F0" }}>
          <div className="text-xs space-y-1" style={{ color: "#64748B" }}>
            <div style={{ color: "#334155", fontWeight: 500 }}>Aaryan Senthilvanan · S.Y.A.L.I.S Labs · 2026</div>
            <div>All data, code, and generated structures are openly available.</div>
          </div>
          <div className="flex gap-2">
            <a href="https://github.com/aesenthilvanan-coder/REFOLD"
               target="_blank" rel="noopener noreferrer"
               className="px-4 py-2 border text-xs font-semibold transition-colors"
               style={{ background: "#EFF6FF", borderColor: "#BFDBFE", color: "#1E3A8A", borderRadius: "2px" }}>
              GitHub →
            </a>
            <Link href="/database"
               className="px-4 py-2 border text-xs font-semibold transition-colors"
               style={{ background: "#F8FAFC", borderColor: "#E2E8F0", color: "#334155", borderRadius: "2px" }}>
              PCD Database →
            </Link>
          </div>
        </div>

      </main>
    </div>
  );
}
