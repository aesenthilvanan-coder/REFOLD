"use client";
import { useState, useEffect, useRef, useCallback } from "react";
import { ATLAS_RAW_URL } from "../lib/atlas";

interface StatsData {
  entries: number;
  avgDrug: number;
  diseases: number;
  total: number;   // total ClinVar variants
}

// Ease-out cubic count-up that animates from `from` to `to`
function useAnimatedNumber(to: number, from = 0, duration = 1400) {
  const [display, setDisplay] = useState(from);
  const frameRef = useRef<number>(0);
  const prevTo = useRef(from);

  useEffect(() => {
    const startVal = prevTo.current;
    const startTime = performance.now();

    function tick(now: number) {
      const t = Math.min((now - startTime) / duration, 1);
      const eased = 1 - Math.pow(1 - t, 3); // ease-out cubic
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

function AnimatedBar({ pct, ready }: { pct: number; ready: boolean }) {
  const [width, setWidth] = useState(0);

  useEffect(() => {
    if (!ready) return;
    // Small delay so the CSS transition fires after mount
    const t = setTimeout(() => setWidth(Math.max(pct, 0.08)), 80);
    return () => clearTimeout(t);
  }, [pct, ready]);

  return (
    <div className="score-bar-track">
      <div
        className="score-bar-fill"
        style={{
          width: `${width}%`,
          transition: "width 1.6s cubic-bezier(0.16, 1, 0.3, 1)",
        }}
      />
    </div>
  );
}

export function LiveStats({ initialEntries, initialAvgDrug, initialDiseases, initialTotal }:
  { initialEntries: number; initialAvgDrug: number; initialDiseases: number; initialTotal: number }) {

  const [stats, setStats] = useState<StatsData>({
    entries: initialEntries,
    avgDrug: initialAvgDrug,
    diseases: initialDiseases,
    total: initialTotal,
  });
  const [ready, setReady] = useState(false);

  // On mount, mark ready so animations fire from 0 → initial value
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
      .catch(() => {/* keep showing last known values */});
  }, []);

  // Poll every 30 seconds for new entries
  useEffect(() => {
    const id = setInterval(poll, 30_000);
    return () => clearInterval(id);
  }, [poll]);

  // Animated display values (start from 0 on mount)
  const dispEntries  = useAnimatedNumber(ready ? stats.entries  : 0, 0, 1400);
  const dispAvgDrug  = useAnimatedNumber(ready ? stats.avgDrug  : 0, 0, 1200);
  const dispDiseases = useAnimatedNumber(ready ? stats.diseases : 0, 0, 1300);

  const pct = (stats.entries / stats.total) * 100;

  return (
    <>
      {/* Stat cards */}
      <div className="grid grid-cols-3 gap-3">
        {[
          { v: Math.round(dispEntries).toLocaleString(), l: "Entries" },
          { v: dispAvgDrug.toFixed(3), l: "Avg Drug Score" },
          { v: Math.round(dispDiseases).toLocaleString(), l: "Diseases" },
        ].map(s => (
          <div key={s.l} className="card p-4 text-center">
            <div className="text-2xl font-black font-mono tabular-nums" style={{ color: "var(--accent)" }}>
              {s.v}
            </div>
            <div className="data-label mt-1">{s.l}</div>
          </div>
        ))}
      </div>

      {/* Proteome progress bar */}
      <div className="border border-[var(--border)] rounded-lg px-5 py-3 flex items-center gap-4"
           style={{ background: "var(--bg-2)" }}>
        <div className="flex-1">
          <div className="flex justify-between items-baseline mb-1.5">
            <span className="section-label">Proteome Coverage</span>
            <span className="text-xs font-mono tabular-nums" style={{ color: "var(--text-secondary)" }}>
              {Math.round(dispEntries).toLocaleString()} / {stats.total.toLocaleString()} variants ({pct.toFixed(4)}%)
            </span>
          </div>
          <AnimatedBar pct={pct} ready={ready} />
        </div>
        <div className="flex items-center gap-1.5 shrink-0">
          <div className="w-1.5 h-1.5 rounded-full bg-[var(--green)] pulse-dot" />
          <span className="text-[10px] font-mono" style={{ color: "var(--green)" }}>LIVE</span>
        </div>
      </div>
    </>
  );
}
