"use client";
import { useState, useEffect, useRef, useCallback } from "react";
import { ATLAS_RAW_URL } from "../lib/atlas";

interface StatsData {
  entries: number;
  avgDrug: number;
  diseases: number;
  total: number;
}

function useAnimatedNumber(to: number, from = 0, duration = 1000) {
  const [display, setDisplay] = useState(from);
  const frameRef = useRef<number>(0);
  const prevTo = useRef(from);

  useEffect(() => {
    const startVal = prevTo.current;
    const startTime = performance.now();
    function tick(now: number) {
      const t = Math.min((now - startTime) / duration, 1);
      const eased = 1 - Math.pow(1 - t, 3);
      setDisplay(startVal + (to - startVal) * eased);
      if (t < 1) frameRef.current = requestAnimationFrame(tick);
      else prevTo.current = to;
    }
    cancelAnimationFrame(frameRef.current);
    frameRef.current = requestAnimationFrame(tick);
    return () => cancelAnimationFrame(frameRef.current);
  }, [to, duration]);

  return display;
}

export function LiveStats({ initialEntries, initialAvgDrug, initialDiseases, initialTotal }:
  { initialEntries: number; initialAvgDrug: number; initialDiseases: number; initialTotal: number }) {

  const [stats, setStats] = useState<StatsData>({
    entries: initialEntries, avgDrug: initialAvgDrug,
    diseases: initialDiseases, total: initialTotal,
  });
  const [ready, setReady] = useState(false);
  useEffect(() => { setReady(true); }, []);

  const poll = useCallback(() => {
    fetch(ATLAS_RAW_URL, { cache: "no-store" })
      .then(r => r.json())
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      .then((data: any) => {
        const entries: number = data.total_entries ?? 0;
        const avgDrug: number = entries > 0
          ? data.entries.reduce((s: number, e: any) => s + (e.pocket?.fpocket_druggability ?? 0), 0) / entries
          : 0;
        const diseases: number = new Set(data.entries.map((e: any) => e.metadata?.disease)).size;
        const total: number = data.proteome_targets?.total_clinvar_pathogenic_missense ?? 178597;
        setStats({ entries, avgDrug, diseases, total });
      })
      .catch(() => {});
  }, []);

  useEffect(() => {
    const id = setInterval(poll, 30_000);
    return () => clearInterval(id);
  }, [poll]);

  const dispEntries  = useAnimatedNumber(ready ? stats.entries  : 0, 0, 1000);
  const dispAvgDrug  = useAnimatedNumber(ready ? stats.avgDrug  : 0, 0, 900);
  const dispDiseases = useAnimatedNumber(ready ? stats.diseases : 0, 0, 950);
  const pct = (stats.entries / stats.total) * 100;

  return (
    <div className="space-y-2">
      {/* Compact stat table */}
      <div className="border" style={{ borderColor: "var(--border)", overflow: "hidden" }}>
        <table className="w-full border-collapse text-xs">
          <thead>
            <tr style={{ background: "var(--bg-2)", borderBottom: "1px solid var(--border)" }}>
              <th className="px-3 py-1.5 text-left font-semibold uppercase tracking-wider"
                  style={{ color: "var(--text-muted)", fontSize: 10 }}>Metric</th>
              <th className="px-3 py-1.5 text-right font-semibold uppercase tracking-wider"
                  style={{ color: "var(--text-muted)", fontSize: 10 }}>Value</th>
            </tr>
          </thead>
          <tbody>
            {[
              { label: "Total entries", value: Math.round(dispEntries).toLocaleString(), color: "var(--accent)" },
              { label: "Avg. fpocket drug score", value: dispAvgDrug.toFixed(3), color: "var(--accent)" },
              { label: "Diseases represented", value: Math.round(dispDiseases).toLocaleString(), color: "var(--text-primary)" },
              { label: "ClinVar pathogenic missense (target)", value: stats.total.toLocaleString(), color: "var(--text-muted)" },
            ].map((row, i) => (
              <tr key={row.label}
                  style={{ background: i % 2 === 0 ? "var(--bg)" : "var(--bg-1)", borderBottom: "1px solid var(--border)" }}>
                <td className="px-3 py-1.5 text-xs" style={{ color: "var(--text-secondary)" }}>{row.label}</td>
                <td className="px-3 py-1.5 text-right font-mono font-semibold tabular-nums"
                    style={{ color: row.color }}>{row.value}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      {/* Proteome coverage bar */}
      <div>
        <div className="flex justify-between items-baseline mb-1">
          <span className="section-label">Proteome coverage</span>
          <div className="flex items-center gap-2">
            <span className="text-[10px] font-mono tabular-nums" style={{ color: "var(--text-muted)" }}>
              {pct.toFixed(4)}%
            </span>
            <div className="flex items-center gap-1">
              <div className="w-1 h-1 rounded-full pulse-dot" style={{ background: "var(--green)" }} />
              <span className="text-[10px] font-mono" style={{ color: "var(--green)", fontWeight: 700 }}>LIVE</span>
            </div>
          </div>
        </div>
        <div className="score-bar-track">
          <div className="score-bar-fill" style={{ width: `${Math.max(pct, 0.05)}%`,
            transition: "width 1.4s cubic-bezier(0.16, 1, 0.3, 1)" }} />
        </div>
      </div>
    </div>
  );
}
