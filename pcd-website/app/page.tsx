import Link from "next/link";
import Image from "next/image";
import { PCDAtlas } from "./types";
import { fetchAtlas } from "./lib/atlas";
import { GnnDiagram } from "./components/GnnDiagram";
import { LiveStats } from "./components/LiveStats";

export const dynamic = "force-dynamic";

// ── Nav ─────────────────────────────────────────────────────────────────────
function NavBar() {
  return (
    <nav className="fixed top-0 left-0 right-0 z-50 border-b border-[var(--border)] backdrop-blur-md"
         style={{ background: "rgba(7,9,15,0.92)" }}>
      <div className="max-w-7xl mx-auto px-6 flex items-center justify-between h-14">
        <div className="flex items-center gap-3">
          <div className="w-7 h-7 rounded-md border border-[var(--accent-border)] flex items-center justify-center"
               style={{ background: "var(--accent-dim)" }}>
            <span className="text-[var(--accent)] font-bold text-xs font-mono">P</span>
          </div>
          <span className="font-semibold text-[var(--text-primary)] text-sm tracking-wide">PCD</span>
          <span className="text-[var(--text-muted)] text-xs hidden sm:inline">· S.Y.A.L.I.S Labs</span>
        </div>
        <div className="flex items-center gap-6 text-sm">
          {["science","architecture","database","about"].map(s => (
            <a key={s} href={`#${s}`}
               className="text-[var(--text-secondary)] hover:text-[var(--text-primary)] transition-colors capitalize text-xs tracking-wide hidden md:block">
              {s}
            </a>
          ))}
          <Link href="/database"
            className="px-3 py-1.5 rounded-md text-xs font-semibold border transition-colors"
            style={{ borderColor: "var(--accent-border)", color: "var(--accent)", background: "var(--accent-dim)" }}>
            Explore →
          </Link>
        </div>
      </div>
    </nav>
  );
}

// ── Hero ─────────────────────────────────────────────────────────────────────
function HeroSection({ data }: { data: PCDAtlas }) {
  const avgDrug = data.entries.reduce((s, e) => s + (e.pocket?.fpocket_druggability ?? 0), 0) / data.entries.length;
  const diseases = new Set(data.entries.map(e => e.metadata.disease)).size;
  const total = data.proteome_targets?.total_clinvar_pathogenic_missense ?? 178597;

  return (
    <section className="relative min-h-screen flex flex-col justify-center pt-14 overflow-hidden">
      <div className="absolute inset-0 grid-bg pointer-events-none" />
      <div className="relative z-10 max-w-7xl mx-auto px-6 w-full py-20">
        <div className="grid lg:grid-cols-2 gap-16 items-center">

          <div className="space-y-8">
            <div>
              <div className="flex items-center gap-2 mb-5">
                <span className="tag tag-green">World First</span>
                <span className="tag tag-cyan">v1.0 — Live</span>
              </div>
              <h1 className="text-5xl font-black tracking-tight leading-[1.1]" style={{ color: "var(--text-primary)" }}>
                Pharmacological<br />
                <span style={{ color: "var(--accent)" }}>Chaperone</span><br />
                Database
              </h1>
              <p className="mt-5 text-base leading-7" style={{ color: "var(--text-secondary)" }}>
                The world's first proteome-scale computational database of de novo
                pharmacological chaperones for pathogenic protein misfolding variants.
                Built by <span style={{ color: "var(--text-primary)" }}>Aaryan Senthilvanan</span> using REFOLD.
              </p>
            </div>

            <LiveStats
              initialEntries={data.total_entries}
              initialAvgDrug={avgDrug}
              initialDiseases={diseases}
              initialTotal={total}
            />

            <div className="flex gap-3">
              <Link href="/database"
                className="px-5 py-2.5 rounded-lg text-sm font-semibold border transition-colors"
                style={{ background: "var(--accent-dim)", borderColor: "var(--accent-border)", color: "var(--accent)" }}>
                Browse Database →
              </Link>
              <a href="#science"
                className="px-5 py-2.5 rounded-lg text-sm font-semibold border transition-colors"
                style={{ borderColor: "var(--border)", color: "var(--text-secondary)" }}>
                How It Works
              </a>
            </div>
          </div>

          <div className="flex flex-col gap-6 items-center">
            <div className="relative w-full max-w-md">
              <div className="rounded-xl overflow-hidden border border-[var(--border)]" style={{ background: "#000" }}>
                <Image
                  src="/images/cftr_molecule.png"
                  alt="3D molecular rendering of a pharmacological chaperone"
                  width={600}
                  height={450}
                  className="w-full h-auto"
                  priority
                />
              </div>
              <div className="absolute -bottom-2 -right-2 tag tag-cyan text-[9px]">
                REFOLD-PC-CFTR-001
              </div>
            </div>
          </div>
        </div>
      </div>
    </section>
  );
}

// ── Science section ────────────────────────────────────────────────────────
function ScienceSection() {
  return (
    <section id="science" className="py-20 border-t border-[var(--border)]">
      <div className="max-w-7xl mx-auto px-6 space-y-16">

        <div className="max-w-2xl">
          <div className="section-label mb-3">The Problem</div>
          <h2 className="text-3xl font-black" style={{ color: "var(--text-primary)" }}>
            Protein Misfolding & The Chaperone Gap
          </h2>
          <p className="mt-4 leading-7" style={{ color: "var(--text-secondary)" }}>
            40–60% of pathogenic missense variants cause disease not by disrupting a protein's
            active site, but by causing it to misfold — the endoplasmic reticulum quality
            control machinery flags it as defective and routes it to proteasomal degradation
            before it can reach its functional destination.
          </p>
          <p className="mt-3 leading-7" style={{ color: "var(--text-secondary)" }}>
            The key insight: <span style={{ color: "var(--text-primary)", fontWeight: 600 }}>
            these proteins work. They just never arrive.
            </span> A small molecule that binds the misfolded intermediate and shifts the folding
            equilibrium back toward the native state is a pharmacological chaperone — and no
            computational framework for designing them has ever been built. REFOLD is that framework.
          </p>
        </div>

        {/* Chaperone Gap + Pipeline diagram (Image #2) */}
        <div className="space-y-3">
          <div className="section-label">REFOLD Architecture — Systemic Overview</div>
          <div className="rounded-xl overflow-hidden border border-[var(--border)]">
            <Image
              src="/images/arch_chaperone_gap.png"
              alt="The Pharmacological Chaperone Gap and REFOLD Pipeline systemic overview"
              width={1400}
              height={800}
              className="w-full h-auto"
            />
          </div>
          <p className="text-xs" style={{ color: "var(--text-muted)" }}>
            Top: Native vs. mutant folding energy landscape. Bottom: REFOLD Stage 1→2→3 pipeline
            producing the proteome-wide chaperone map.
          </p>
        </div>

        {/* Three-stage pipeline cards */}
        <div className="grid md:grid-cols-3 gap-4">
          {[
            {
              stage: "Stage 1", name: "Rescue Amenability Classification",
              label: "GNN-LM Fusion Classifier", color: "var(--accent)",
              points: [
                "AlphaFold2 structure + ESM-2 embeddings (1280-d)",
                "FoldX ΔΔG + evolutionary conservation PSSM",
                "Atom-edge GNN message passing over protein graph",
                "AUC 0.88 (Gaucher), 0.79 (Cystic Fibrosis)",
              ],
            },
            {
              stage: "Stage 2", name: "Conformational Sampling & Pocket ID",
              label: "ANM + fpocket", color: "var(--violet)",
              points: [
                "Anisotropic Network Model over AlphaFold structure",
                "25–50 conformations sampled along non-trivial modes",
                "fpocket α-sphere druggability scoring per conformation",
                "Cryptic pocket isolation: absent in WT, drug > 0.75",
              ],
            },
            {
              stage: "Stage 3", name: "De Novo Drug Design",
              label: "SE(3) Equivariant Diffusion", color: "var(--green)",
              points: [
                "E_ij pocket geometry matrix as conditioning input",
                "SE(3)-equivariant E(3)-equivariant denoising network",
                "T=100 → T=1 diffusion on atom coordinates",
                "Composite score: 0.45×affinity + 0.30×QED + 0.15×SA + 0.10×logP",
              ],
            },
          ].map(s => (
            <div key={s.stage} className="card p-5 space-y-4">
              <span className="tag text-[10px]"
                    style={{ background: `${s.color}15`, color: s.color, borderColor: `${s.color}30` }}>
                {s.stage}
              </span>
              <div>
                <div className="font-bold text-sm" style={{ color: "var(--text-primary)" }}>{s.name}</div>
                <div className="text-xs mt-0.5 font-mono" style={{ color: s.color }}>{s.label}</div>
              </div>
              <ul className="space-y-1.5">
                {s.points.map((p, i) => (
                  <li key={i} className="flex gap-2 text-xs leading-5" style={{ color: "var(--text-secondary)" }}>
                    <span style={{ color: s.color, flexShrink: 0 }}>·</span>{p}
                  </li>
                ))}
              </ul>
            </div>
          ))}
        </div>

        {/* GNN-LM + SE(3) diagram (Image #1 / #3) */}
        <div id="architecture" className="space-y-3">
          <div className="section-label">ML Architecture — Technical Detail</div>
          <GnnDiagram />
          <p className="text-xs" style={{ color: "var(--text-muted)" }}>
            Part A: GNN feature extractor with atom-edge message passing + Transformer LM fusion → rescuability logits.
            Part B: SE(3)-equivariant denoising network for structure-conditioned molecular generation (T=T→T=1).
          </p>
        </div>

        {/* Pipeline overview diagram */}
        <div className="space-y-3">
          <div className="section-label">Pipeline — End-to-End Flow</div>
          <div className="rounded-xl overflow-hidden border border-[var(--border)]">
            <Image
              src="/images/arch_pipeline_overview.png"
              alt="REFOLD pipeline end-to-end overview"
              width={1400}
              height={700}
              className="w-full h-auto"
            />
          </div>
        </div>

      </div>
    </section>
  );
}

// ── Database preview ────────────────────────────────────────────────────────
function DatabasePreview({ data }: { data: PCDAtlas }) {
  const preview = data.entries.slice(0, 3);
  return (
    <section id="database" className="py-20 border-t border-[var(--border)]">
      <div className="max-w-7xl mx-auto px-6 space-y-8">
        <div className="flex items-end justify-between">
          <div>
            <div className="section-label mb-2">Live Atlas</div>
            <h2 className="text-3xl font-black" style={{ color: "var(--text-primary)" }}>Chaperone Database</h2>
            <p className="mt-2 text-sm" style={{ color: "var(--text-secondary)" }}>
              {data.total_entries} entries · continuously growing · powered by REFOLD
            </p>
          </div>
          <Link href="/database"
            className="px-4 py-2 rounded-lg text-sm font-semibold border transition-colors"
            style={{ borderColor: "var(--border)", color: "var(--text-secondary)" }}>
            View all →
          </Link>
        </div>

        <div className="space-y-2">
          {preview.map(entry => (
            <Link key={entry.entry_id} href={`/database/${entry.entry_id}`} className="block">
              <div className="card-hover rounded-lg p-4 grid grid-cols-12 gap-4 items-center">
                <div className="col-span-3 space-y-1">
                  <div className="flex items-center gap-2">
                    <span className="tag tag-cyan">{entry.metadata.gene}</span>
                    <span className="tag tag-amber">{entry.metadata.mutation_mature}</span>
                  </div>
                  <div className="text-[10px] font-mono" style={{ color: "var(--text-muted)" }}>{entry.entry_id}</div>
                </div>
                <div className="col-span-4">
                  <div className="text-sm font-medium" style={{ color: "var(--text-primary)" }}>{entry.metadata.disease}</div>
                  <div className="text-xs mt-0.5 line-clamp-1" style={{ color: "var(--text-muted)" }}>{entry.metadata.mechanism}</div>
                </div>
                <div className="col-span-3 space-y-1.5">
                  <div className="flex justify-between text-[10px] mb-0.5">
                    <span style={{ color: "var(--text-muted)" }}>Drug Score</span>
                    <span className="font-mono" style={{ color: "var(--accent)" }}>{entry.pocket.fpocket_druggability.toFixed(3)}</span>
                  </div>
                  <div className="score-bar-track">
                    <div className="score-bar-fill" style={{ width: `${entry.pocket.fpocket_druggability * 100}%` }} />
                  </div>
                </div>
                <div className="col-span-2 flex justify-end">
                  <span className="text-xs" style={{ color: "var(--accent)" }}>→</span>
                </div>
              </div>
            </Link>
          ))}
        </div>

        <div className="text-center">
          <Link href="/database"
            className="inline-block px-6 py-2.5 rounded-lg text-sm font-semibold border transition-colors"
            style={{ borderColor: "var(--accent-border)", color: "var(--accent)", background: "var(--accent-dim)" }}>
            Browse Full Database
          </Link>
        </div>
      </div>
    </section>
  );
}

// ── About ────────────────────────────────────────────────────────────────────
function AboutSection() {
  return (
    <section id="about" className="py-20 border-t border-[var(--border)]">
      <div className="max-w-7xl mx-auto px-6">
        <div className="grid lg:grid-cols-3 gap-12 items-start">

          <div className="space-y-5">
            <div className="rounded-xl overflow-hidden border border-[var(--border)] w-48">
              <Image
                src="/images/aaryan.jpg"
                alt="Aaryan Senthilvanan"
                width={300}
                height={400}
                className="w-full h-auto"
              />
            </div>
            <div>
              <div className="font-bold text-base" style={{ color: "var(--text-primary)" }}>Aaryan Senthilvanan</div>
              <div className="text-sm mt-0.5" style={{ color: "var(--text-muted)" }}>Investigator</div>
              <div className="text-sm mt-0.5 font-mono" style={{ color: "var(--accent)" }}>S.Y.A.L.I.S Labs</div>
            </div>
          </div>

          <div className="lg:col-span-2 space-y-6">
            <div>
              <div className="section-label mb-3">About</div>
              <h2 className="text-3xl font-black" style={{ color: "var(--text-primary)" }}>
                Building the cure for every misfolding disease
              </h2>
            </div>

            <div className="space-y-4 text-sm leading-7" style={{ color: "var(--text-secondary)" }}>
              <p>
                REFOLD is a fully automated, end-to-end computational pipeline that takes any
                human disease-associated missense mutation as input and outputs: (1) a binary
                prediction of whether the mutation causes misfolding rather than functional
                disruption; (2) the 3D structure of the mutant protein in its partially unfolded
                state with transient binding pockets identified; (3) de novo designed small molecule
                candidates that bind those pockets and thermodynamically stabilize the correct fold.
              </p>
              <p>
                The PCD is the output of running REFOLD continuously across the entire human
                proteome — every pathogenic missense variant in ClinVar, processed one by one,
                with each result injected into this living database automatically.
              </p>
              <blockquote className="border-l-2 pl-4 py-1" style={{ borderColor: "var(--accent)" }}>
                <p className="italic" style={{ color: "var(--text-primary)" }}>
                  "Three FDA-approved pharmacological chaperones exist. Each is a multi-billion
                  dollar drug discovered by brute-force screening of hundreds of thousands of
                  compounds against a single protein. REFOLD designs them from first principles,
                  for every protein, automatically."
                </p>
              </blockquote>
            </div>

            <div className="grid grid-cols-3 gap-3">
              {[
                { v: "20,430", l: "Canonical human proteins" },
                { v: "178,597", l: "ClinVar pathogenic missense variants" },
                { v: "5,992", l: "Genes with known pathogenic variants" },
              ].map(s => (
                <div key={s.l} className="card p-4">
                  <div className="text-xl font-black font-mono" style={{ color: "var(--accent)" }}>{s.v}</div>
                  <div className="data-label mt-1">{s.l}</div>
                </div>
              ))}
            </div>
          </div>
        </div>
      </div>
    </section>
  );
}

// ── Footer ───────────────────────────────────────────────────────────────────
function Footer({ data }: { data: PCDAtlas }) {
  return (
    <footer className="border-t border-[var(--border)] py-8">
      <div className="max-w-7xl mx-auto px-6 flex flex-col sm:flex-row items-center justify-between gap-4">
        <div className="flex items-center gap-3">
          <div className="w-6 h-6 rounded border border-[var(--accent-border)] flex items-center justify-center"
               style={{ background: "var(--accent-dim)" }}>
            <span className="text-[var(--accent)] font-bold text-[10px] font-mono">P</span>
          </div>
          <span className="text-xs" style={{ color: "var(--text-muted)" }}>
            PCD v{data.version} · S.Y.A.L.I.S Labs · Powered by REFOLD
          </span>
        </div>
        <div className="text-xs font-mono" style={{ color: "var(--text-muted)" }}>
          {data.total_entries} entries · {data.build_date}
        </div>
      </div>
    </footer>
  );
}

// ── Page (server component, fetches live data) ───────────────────────────────
export default async function HomePage() {
  const data = await fetchAtlas();

  return (
    <div style={{ background: "var(--bg)", minHeight: "100vh" }}>
      <NavBar />
      <HeroSection data={data} />
      <ScienceSection />
      <DatabasePreview data={data} />
      <AboutSection />
      <Footer data={data} />
    </div>
  );
}
