"use client";
import { useState, useMemo, useEffect } from "react";
import Link from "next/link";
import { useRouter } from "next/navigation";
import { ATLAS_RAW_URL } from "../lib/atlas";
import { PCDAtlas } from "../types";

// eslint-disable-next-line @typescript-eslint/no-explicit-any
function EntryRow({ entry, index }: { entry: any; index: number }) {
  const router = useRouter();
  const meta = entry.metadata ?? {};
  const pkt  = entry.pocket   ?? {};
  const chap = entry.chaperone ?? {};
  const drug = pkt.fpocket_druggability ?? 0;
  const score = chap.composite_score ?? 0;

  return (
    <tr
      onClick={() => router.push(`/database/${entry.entry_id}`)}
      className="cursor-pointer"
      style={{
        background: index % 2 === 0 ? "var(--bg)" : "var(--bg-1)",
        borderBottom: "1px solid var(--border)",
      }}
      onMouseOver={e => ((e.currentTarget as HTMLElement).style.background = "var(--accent-dim)")}
      onMouseOut={e => ((e.currentTarget as HTMLElement).style.background = index % 2 === 0 ? "var(--bg)" : "var(--bg-1)")}
    >
      <td className="px-4 py-2.5 whitespace-nowrap">
        <div className="flex items-center gap-1.5">
          <span className="tag tag-cyan">{meta.gene}</span>
          <span className="tag tag-amber">{meta.mutation_mature}</span>
        </div>
        <div className="font-mono text-[10px] mt-0.5" style={{ color: "var(--text-muted)" }}>{entry.entry_id}</div>
      </td>
      <td className="px-4 py-2.5 max-w-xs">
        <div className="text-xs font-medium" style={{ color: "var(--text-primary)" }}>{meta.disease}</div>
        <div className="text-[11px] mt-0.5 truncate" style={{ color: "var(--text-muted)" }}>{meta.mechanism}</div>
      </td>
      <td className="px-4 py-2.5 whitespace-nowrap">
        <div className="font-mono font-semibold text-xs tabular-nums" style={{ color: "var(--accent)" }}>
          {drug.toFixed(3)}
        </div>
        <div className="mt-1 w-20">
          <div className="score-bar-track">
            <div className="score-bar-fill" style={{ width: `${Math.min(drug * 100, 100)}%` }} />
          </div>
        </div>
      </td>
      <td className="px-4 py-2.5 whitespace-nowrap">
        <div className="font-mono font-semibold text-xs tabular-nums" style={{ color: "var(--violet)" }}>
          {score.toFixed(3)}
        </div>
        <div className="mt-1 w-20">
          <div className="score-bar-track">
            <div className="score-bar-fill" style={{ width: `${Math.min(score * 100, 100)}%`, background: "var(--violet)" }} />
          </div>
        </div>
      </td>
      <td className="px-4 py-2.5 whitespace-nowrap">
        <div className="grid grid-cols-3 gap-3 text-center">
          {[
            { v: (chap.mw  ?? 0).toFixed(0), l: "MW" },
            { v: (chap.logp ?? 0).toFixed(1), l: "logP" },
            { v: (chap.qed  ?? 0).toFixed(2), l: "QED" },
          ].map(p => (
            <div key={p.l}>
              <div className="font-mono font-semibold text-[11px] tabular-nums" style={{ color: "var(--text-primary)" }}>{p.v}</div>
              <div className="data-label">{p.l}</div>
            </div>
          ))}
        </div>
      </td>
      <td className="px-4 py-2.5 text-right">
        <span className="text-xs" style={{ color: "var(--border-mid)" }}>→</span>
      </td>
    </tr>
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

      {/* ── Nav ── */}
      <nav className="fixed top-0 left-0 right-0 z-50 border-b"
           style={{ background: "var(--bg)", borderColor: "var(--border)" }}>
        <div className="max-w-7xl mx-auto px-6 flex items-center gap-3 h-10 text-[11px]">
          <Link href="/" className="font-mono font-bold tracking-widest" style={{ color: "var(--accent)" }}>PCD</Link>
          <span style={{ color: "var(--border)" }}>|</span>
          <span style={{ color: "var(--text-secondary)", fontWeight: 500 }}>Chaperone Database</span>
          <div className="ml-auto flex items-center gap-4">
            <Link href="/abstract" className="nav-link">Abstract</Link>
            <Link href="/" className="nav-link">Home</Link>
            {atlas && (
              <div className="flex items-center gap-1.5">
                <div className="w-1.5 h-1.5 rounded-full pulse-dot" style={{ background: "var(--green)" }} />
                <span className="font-mono text-[10px]" style={{ color: "var(--green)" }}>LIVE</span>
              </div>
            )}
          </div>
        </div>
      </nav>

      <main className="max-w-7xl mx-auto px-6 pb-12" style={{ paddingTop: "52px" }}>

        {!atlas ? (
          <div className="flex items-center justify-center py-32">
            <div className="text-center space-y-2">
              <div className="flex justify-center gap-1">
                {[0,1,2].map(i => (
                  <div key={i} className="w-1.5 h-1.5 rounded-full pulse-dot"
                       style={{ background: "var(--accent)", animationDelay: `${i*0.2}s` }} />
                ))}
              </div>
              <div className="section-label">Loading atlas...</div>
            </div>
          </div>
        ) : (
          <>
            {/* ── Page header ── */}
            <div className="border-b py-5 mb-5" style={{ borderColor: "var(--border)" }}>
              <div className="flex flex-col sm:flex-row gap-4 items-start sm:items-end justify-between">
                <div>
                  <div className="section-label mb-1">Pharmacological Chaperone Database</div>
                  <h1 className="text-lg font-bold" style={{ color: "var(--text-primary)" }}>Entry Browser</h1>
                  <p className="text-xs mt-0.5 font-mono" style={{ color: "var(--text-muted)" }}>
                    {atlas.total_entries.toLocaleString()} entries ·
                    {" "}{new Set(atlas.entries.map(e => e.metadata.disease)).size} diseases ·
                    {" "}avg. drug score {(atlas.entries.reduce((s: number, e) => s + (e.pocket?.fpocket_druggability ?? 0), 0) / atlas.entries.length).toFixed(3)} ·
                    {" "}build {atlas.build_date}
                  </p>
                </div>

                {/* Proteome progress */}
                {atlas.proteome_targets && (
                  <div className="shrink-0 w-full sm:w-72">
                    <div className="flex justify-between items-baseline mb-1">
                      <span className="section-label">Proteome Coverage</span>
                      <span className="text-[10px] font-mono" style={{ color: "var(--text-muted)" }}>
                        {atlas.total_entries.toLocaleString()} / {atlas.proteome_targets.total_clinvar_pathogenic_missense.toLocaleString()}
                      </span>
                    </div>
                    <div className="score-bar-track">
                      <div className="score-bar-fill"
                           style={{ width: `${Math.max((atlas.total_entries / atlas.proteome_targets.total_clinvar_pathogenic_missense) * 100, 0.4)}%` }} />
                    </div>
                  </div>
                )}
              </div>
            </div>

            {/* ── Search & Sort ── */}
            <div className="flex flex-col sm:flex-row gap-2 mb-3">
              <input
                type="text"
                placeholder="Search gene, disease, mutation, SMILES..."
                value={query}
                onChange={e => setQuery(e.target.value)}
                className="flex-1 px-3 py-1.5 text-xs border outline-none"
                style={{
                  background: "var(--bg)",
                  borderColor: "var(--border-mid)",
                  color: "var(--text-primary)",
                  borderRadius: "2px",
                  fontFamily: "inherit",
                }}
                onFocus={e => (e.currentTarget.style.borderColor = "var(--accent)")}
                onBlur={e => (e.currentTarget.style.borderColor = "var(--border-mid)")}
              />
              <div className="flex gap-1 shrink-0">
                <span className="flex items-center text-[10px] px-2" style={{ color: "var(--text-muted)" }}>Sort by:</span>
                {(["fpocket","score","mw"] as const).map(k => (
                  <button key={k} onClick={() => setSortKey(k)}
                    className="px-2.5 py-1.5 text-[11px] font-semibold border"
                    style={{
                      background: sortKey === k ? "var(--accent)" : "var(--bg)",
                      color: sortKey === k ? "#fff" : "var(--text-muted)",
                      borderColor: sortKey === k ? "var(--accent)" : "var(--border-mid)",
                      borderRadius: "2px",
                    }}>
                    {k === "fpocket" ? "Drug Score ↓" : k === "score" ? "Chaperone ↓" : "MW ↑"}
                  </button>
                ))}
              </div>
            </div>

            <div className="text-[11px] mb-2 font-mono" style={{ color: "var(--text-muted)" }}>
              Showing {filtered.length.toLocaleString()} of {atlas.total_entries.toLocaleString()} entries
              {query && <span style={{ color: "var(--text-secondary)" }}> — filtered by &ldquo;{query}&rdquo;</span>}
            </div>

            {/* ── Data Table ── */}
            <div className="border" style={{ borderColor: "var(--border)", overflow: "hidden" }}>
              <table className="w-full border-collapse text-xs">
                <thead>
                  <tr style={{ background: "var(--bg-2)", borderBottom: "2px solid var(--border-mid)" }}>
                    {[
                      { label: "Gene / Variant", w: "20%" },
                      { label: "Disease / Mechanism", w: "28%" },
                      { label: "Drug Score", w: "14%", note: "fpocket" },
                      { label: "Chaperone Score", w: "14%", note: "composite" },
                      { label: "Properties", w: "18%", note: "MW · logP · QED" },
                      { label: "", w: "4%" },
                    ].map(h => (
                      <th key={h.label}
                          className="px-4 py-2 text-left font-semibold"
                          style={{ color: "var(--text-muted)", fontSize: 10, width: h.w }}>
                        <span className="uppercase tracking-wider">{h.label}</span>
                        {"note" in h && h.note && (
                          <span className="font-normal ml-1 normal-case tracking-normal" style={{ color: "var(--border-mid)" }}>
                            ({h.note})
                          </span>
                        )}
                      </th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {filtered.map((entry, i) => (
                    <EntryRow key={entry.entry_id} entry={entry} index={i} />
                  ))}
                  {filtered.length === 0 && (
                    <tr>
                      <td colSpan={6} className="text-center py-16 font-mono text-xs"
                          style={{ color: "var(--text-muted)" }}>
                        No entries match &quot;{query}&quot;
                      </td>
                    </tr>
                  )}
                </tbody>
              </table>
            </div>
          </>
        )}
      </main>
    </div>
  );
}
