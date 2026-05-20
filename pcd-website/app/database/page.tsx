"use client";
import { useState, useMemo, useEffect } from "react";
import Link from "next/link";
import { ATLAS_RAW_URL } from "../lib/atlas";
import { PCDAtlas, PCDEntry } from "../types";

function ScoreBar({ value, color = "var(--accent)" }: { value: number; color?: string }) {
  return (
    <div className="flex items-center gap-2">
      <div className="score-bar-track flex-1">
        <div className="score-bar-fill" style={{ width: `${Math.min(value * 100, 100)}%`, background: color }} />
      </div>
      <span className="text-xs font-mono font-bold w-12 text-right" style={{ color }}>{value.toFixed(3)}</span>
    </div>
  );
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any
function EntryRow({ entry }: { entry: any }) {
  const meta = entry.metadata ?? {};
  const pkt  = entry.pocket   ?? {};
  const chap = entry.chaperone ?? {};
  return (
    <Link href={`/database/${entry.entry_id}`} className="block group">
      <div className="card-hover rounded-lg p-4 grid grid-cols-12 gap-4 items-center">
        <div className="col-span-12 sm:col-span-3 space-y-1.5">
          <div className="flex items-center gap-2 flex-wrap">
            <span className="tag tag-cyan">{meta.gene}</span>
            <span className="tag tag-amber">{meta.mutation_mature}</span>
          </div>
          <div className="text-[10px] font-mono" style={{ color: "var(--text-muted)" }}>{entry.entry_id}</div>
        </div>
        <div className="col-span-12 sm:col-span-3">
          <div className="text-sm font-medium" style={{ color: "var(--text-primary)" }}>{meta.disease}</div>
          <div className="text-xs mt-0.5 line-clamp-1" style={{ color: "var(--text-muted)" }}>{meta.mechanism}</div>
        </div>
        <div className="col-span-12 sm:col-span-3 space-y-2">
          <div>
            <div className="data-label mb-1">fpocket Drug.</div>
            <ScoreBar value={pkt.fpocket_druggability ?? 0} />
          </div>
          <div>
            <div className="data-label mb-1">Chaperone Score</div>
            <ScoreBar value={chap.composite_score ?? 0} color="var(--violet)" />
          </div>
        </div>
        <div className="col-span-12 sm:col-span-2 grid grid-cols-3 gap-2 text-center">
          {[
            { v: (chap.mw ?? 0).toFixed(0), l: "MW" },
            { v: (chap.logp ?? 0).toFixed(1), l: "logP" },
            { v: (chap.qed ?? 0).toFixed(2), l: "QED" },
          ].map(p => (
            <div key={p.l}>
              <div className="text-xs font-bold font-mono" style={{ color: "var(--text-primary)" }}>{p.v}</div>
              <div className="data-label">{p.l}</div>
            </div>
          ))}
        </div>
        <div className="col-span-12 sm:col-span-1 flex justify-end">
          <span className="text-xs group-hover:translate-x-0.5 transition-transform"
                style={{ color: "var(--accent)" }}>→</span>
        </div>
      </div>
    </Link>
  );
}

export default function DatabasePage() {
  const [atlas, setAtlas] = useState<PCDAtlas | null>(null);
  const [query, setQuery] = useState("");
  const [sortKey, setSortKey] = useState<"fpocket" | "score" | "mw">("fpocket");

  useEffect(() => {
    let mounted = true;
    const load = () =>
      fetch(ATLAS_RAW_URL, { cache: "no-store" })
        .then(r => r.json())
        .then(data => { if (mounted) setAtlas(data); })
        .catch(console.error);

    load();
    // Poll every 30 s — shows new entries without a manual refresh
    const timer = setInterval(load, 30_000);
    return () => { mounted = false; clearInterval(timer); };
  }, []);

  const filtered = useMemo(() => {
    if (!atlas) return [];
    const q = query.toLowerCase();
    return atlas.entries
      .filter(e =>
        !q ||
        e.metadata.gene.toLowerCase().includes(q) ||
        e.metadata.disease.toLowerCase().includes(q) ||
        e.metadata.mutation_mature.toLowerCase().includes(q) ||
        e.entry_id.toLowerCase().includes(q) ||
        e.chaperone.smiles.toLowerCase().includes(q)
      )
      .sort((a, b) => {
        if (sortKey === "fpocket") return b.pocket.fpocket_druggability - a.pocket.fpocket_druggability;
        if (sortKey === "score")   return b.chaperone.composite_score - a.chaperone.composite_score;
        if (sortKey === "mw")      return a.chaperone.mw - b.chaperone.mw;
        return 0;
      });
  }, [atlas, query, sortKey]);

  return (
    <div style={{ background: "var(--bg)", minHeight: "100vh" }}>
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
          <span className="font-semibold" style={{ color: "var(--text-primary)" }}>Database</span>
          <div className="ml-auto flex items-center gap-4">
            <Link href="/abstract" className="text-xs transition-colors"
                  style={{ color: "var(--text-muted)" }}>Abstract</Link>
            {atlas && (
              <div className="flex items-center gap-1.5">
                <div className="w-1.5 h-1.5 rounded-full bg-[var(--green)] pulse-dot" />
                <span className="text-[10px] font-mono" style={{ color: "var(--green)" }}>LIVE</span>
              </div>
            )}
          </div>
        </div>
      </nav>

      <main className="max-w-7xl mx-auto px-6 pt-20 pb-16 space-y-6">

        {!atlas ? (
          <div className="flex items-center justify-center py-32">
            <div className="text-center space-y-3">
              <div className="flex justify-center gap-1">
                {[0,1,2].map(i => (
                  <div key={i} className="w-1.5 h-1.5 rounded-full bg-[var(--accent)]"
                       style={{ animation: `pulse-dot 1.2s ease-in-out ${i*0.2}s infinite` }} />
                ))}
              </div>
              <div className="section-label">Loading live atlas...</div>
            </div>
          </div>
        ) : (
          <>
            <div className="grid sm:grid-cols-2 gap-6 items-end">
              <div>
                <div className="section-label mb-2">Live Atlas</div>
                <h1 className="text-3xl font-black" style={{ color: "var(--text-primary)" }}>Chaperone Database</h1>
                <p className="text-sm mt-2" style={{ color: "var(--text-secondary)" }}>
                  {atlas.total_entries} entries · Powered by REFOLD · S.Y.A.L.I.S Labs
                </p>
              </div>
              <div className="grid grid-cols-3 gap-3 text-center">
                {[
                  { label: "Total Entries", value: atlas.total_entries },
                  { label: "Avg Drug Score",
                    // eslint-disable-next-line @typescript-eslint/no-explicit-any
                    value: (atlas.entries.reduce((s: number, e: any) => s + (e.pocket?.fpocket_druggability ?? 0), 0) / atlas.entries.length).toFixed(3) },
                  { label: "Diseases",
                    value: new Set(atlas.entries.map(e => e.metadata.disease)).size },
                ].map(s => (
                  <div key={s.label} className="card rounded-lg p-3">
                    <div className="text-xl font-black font-mono" style={{ color: "var(--accent)" }}>{s.value}</div>
                    <div className="data-label mt-1">{s.label}</div>
                  </div>
                ))}
              </div>
            </div>

            {/* Proteome progress */}
            {atlas.proteome_targets && (
              <div className="border border-[var(--border)] rounded-lg px-5 py-3 flex items-center gap-4"
                   style={{ background: "var(--bg-2)" }}>
                <div className="flex-1">
                  <div className="flex justify-between items-baseline mb-1.5">
                    <span className="section-label">Proteome Coverage</span>
                    <span className="text-xs font-mono" style={{ color: "var(--text-secondary)" }}>
                      {atlas.total_entries.toLocaleString()} / {atlas.proteome_targets.total_clinvar_pathogenic_missense.toLocaleString()} variants
                    </span>
                  </div>
                  <div className="score-bar-track">
                    <div className="score-bar-fill" style={{ width: `${Math.max((atlas.total_entries / atlas.proteome_targets.total_clinvar_pathogenic_missense) * 100, 0.5)}%` }} />
                  </div>
                </div>
                <div className="flex items-center gap-1.5 shrink-0">
                  <div className="w-1.5 h-1.5 rounded-full bg-[var(--green)] pulse-dot" />
                  <span className="text-[10px] font-mono" style={{ color: "var(--green)" }}>LIVE</span>
                </div>
              </div>
            )}

            <div className="flex flex-col sm:flex-row gap-2">
              <input
                type="text"
                placeholder="Search gene, disease, mutation, SMILES..."
                value={query}
                onChange={e => setQuery(e.target.value)}
                className="flex-1 px-4 py-2.5 rounded-lg text-sm border outline-none transition-colors"
                style={{ background: "var(--bg-2)", borderColor: "var(--border)", color: "var(--text-primary)" }}
                onFocus={e => (e.currentTarget.style.borderColor = "var(--accent-border)")}
                onBlur={e => (e.currentTarget.style.borderColor = "var(--border)")}
              />
              <div className="flex gap-1.5">
                {(["fpocket","score","mw"] as const).map(k => (
                  <button key={k} onClick={() => setSortKey(k)}
                    className="px-4 py-2.5 rounded-lg text-xs font-semibold border transition-all"
                    style={{
                      background: sortKey===k ? "var(--accent-dim)" : "var(--bg-2)",
                      color: sortKey===k ? "var(--accent)" : "var(--text-secondary)",
                      borderColor: sortKey===k ? "var(--accent-border)" : "var(--border)",
                    }}>
                    {k==="fpocket" ? "Drug Score ↓" : k==="score" ? "Chaperone ↓" : "MW ↑"}
                  </button>
                ))}
              </div>
            </div>

            <div className="text-xs" style={{ color: "var(--text-muted)" }}>
              Showing {filtered.length} of {atlas.total_entries} entries
            </div>

            <div className="space-y-1.5">
              {filtered.map(entry => <EntryRow key={entry.entry_id} entry={entry} />)}
              {filtered.length === 0 && (
                <div className="text-center py-16" style={{ color: "var(--text-muted)" }}>
                  No entries match &quot;{query}&quot;
                </div>
              )}
            </div>
          </>
        )}
      </main>
    </div>
  );
}
